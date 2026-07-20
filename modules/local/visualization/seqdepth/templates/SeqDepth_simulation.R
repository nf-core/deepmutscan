#!/usr/bin/env Rscript

# Sequencing-depth rarefaction: what fraction of the possible library variants stay detectable
# (>= `threshold` reads) as the sequencing depth is subsampled from full down to zero.
#
# This is a closed-form rarefaction, not a stochastic simulation. Subsampling reads without
# replacement makes each variant's surviving read count hypergeometric, so the expected number of
# variants still at or above the threshold at depth n is a sum of hypergeometric tail probabilities:
#
#   E[detected | n reads] = sum_i  P(Hypergeom(N, c_i, n) >= threshold)
#                         = sum_i  phyper(threshold - 1, c_i, N - c_i, n, lower.tail = FALSE)
#
# where N is the total read count and c_i the reads on variant i. It is exact (no Monte-Carlo noise),
# needs no random sampling, and is orders of magnitude faster than draw-one-read-at-a-time loops - so
# the step count can be fine without the runtime blowing up.

library(ggplot2)

SeqDepth_simulation_plot <- function(prefiltered_gatk_path, possible_mutations_path,
                                     output_file_path, curve_csv_path,
                                     n_steps = 40, threshold = 3) {

  data <- read.csv(prefiltered_gatk_path)
  counts <- as.numeric(data\$counts)
  counts <- counts[!is.na(counts) & counts > 0]           # only observed variants carry reads
  possible_mutations <- read.csv(possible_mutations_path)
  total_possible <- nrow(possible_mutations)               # denominator: the whole designed library

  N <- sum(counts)
  if (N <= 0 || total_possible <= 0) {
    stop("SeqDepth: no counts or no possible mutations to rarefy")
  }

  # Evenly spaced depths from 0 to the full library, as integer read totals.
  depths <- unique(round(seq(0, N, length.out = n_steps + 1)))

  # Expected detected variants at each depth. phyper is vectorised over the per-variant counts, so each
  # depth is one vectorised call; `lower.tail = FALSE` gives P(X > threshold-1) = P(X >= threshold).
  expected_detected <- vapply(depths, function(n) {
    if (n <= 0) return(0)
    sum(phyper(threshold - 1, counts, N - counts, n, lower.tail = FALSE))
  }, numeric(1))

  curve <- data.frame(
    depth_fold   = round(depths / N, 4),                       # x: fraction of full sequencing depth
    variants_pct = round(expected_detected / total_possible * 100, 4)  # y: % of the possible library
  )
  write.csv(curve, curve_csv_path, row.names = FALSE)

  plot <- ggplot(curve, aes(x = depth_fold, y = variants_pct)) +
    geom_line(color = "black", linewidth = 0.4) +
    geom_hline(yintercept = 100, linetype = "dotted", color = "black") +
    scale_y_continuous(
      labels = scales::percent_format(scale = 1),
      limits = c(0, 100),
      breaks = seq(0, 100, by = 5)
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, length.out = 20),
      labels = scales::number_format(accuracy = 0.01)
    ) +
    labs(x = "Fraction of sequencing depth", y = "Variants detected (% of possible library)") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid.major = element_line(linewidth = 0.2, linetype = "solid", color = "grey80"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(output_file_path, plot = plot, device = "pdf", width = 8, height = 6)
}


#####
# run function
#####
SeqDepth_simulation_plot(
  prefiltered_gatk_path   = "$variantCounts_filtered_by_library",
  possible_mutations_path = "$possible_mutations",
  output_file_path        = "SeqDepth.pdf",
  curve_csv_path          = "seqdepth_curve.csv",
  n_steps                 = 40,
  threshold               = $min_counts
  )

#####
# create versions.yml
#####
r_version <- strsplit(version[['version.string']], ' ')[[1]][3]
ggplot2_version <- as.character(packageVersion("ggplot2"))

if (is.null(r_version)) r_version <- "unknown"
if (length(ggplot2_version) == 0) ggplot2_version <- "unknown"

f <- file("versions.yml", "w")
writeLines(
  c(
    '"${task.process}":',
    paste('    r-base:', r_version),
    paste('    r-ggplot2:', ggplot2_version)
  ),
  f
)
close(f)
