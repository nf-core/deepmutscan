#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(scales)

# input: completed library table (all programmed variants; counts_per_cov 0 for variants with no
#        counts, observed = FALSE for those never observed)
# output: two LogDiff plots (counts per coverage, sorted, log scale): a line plot, and a dot plot
#         coloured by the number of changed bases per codon
#
# Variants with no counts cannot be drawn on a log scale. They are placed at a floor derived from
# the data - half the lowest non-zero counts per coverage in the library, i.e. below every observed
# variant - and shown in grey. LogDiff is log10 of the 90th over the 10th percentile of all
# programmed variants. When more than 10 % of the library has no counts, the 10th percentile would
# be the floor, an artefact of the plot, so LogDiff is not reported; the share of variants without
# counts, printed on every plot, is then the meaningful measure.

logdiff_plot <- function(completed_path) {

  completed <- read.table(completed_path, sep = ",", fill = NA, header = TRUE)
  cpc <- completed\$counts_per_cov
  cpc[is.na(cpc)] <- 0
  if (!"observed" %in% names(completed)) {
    completed\$observed <- cpc > 0
  }

  n_total    <- length(cpc)
  n_never    <- sum(!completed\$observed)
  n_zeroed   <- sum(completed\$observed & cpc == 0)    # observed, but corrected to 0
  n_no_count <- n_never + n_zeroed

  nonzero <- cpc[cpc > 0]
  if (length(nonzero) == 0) {
    # Nothing to rank: write placeholder plots that say so, instead of failing the run
    empty <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No programmed variant has counts in this library.") +
      theme_void()
    ggsave(filename = "logdiff_plot.pdf", plot = empty, device = "pdf", width = 10, height = 6)
    ggsave(filename = "logdiff_varying_bases.pdf", plot = empty, device = "pdf", width = 10, height = 6)
    return(invisible(NULL))
  }
  floor_value <- min(nonzero) / 2

  sorted_counts <- completed %>%
    mutate(counts_per_cov_plot = ifelse(cpc > 0, cpc, floor_value),
           category = ifelse(cpc > 0, as.character(varying_bases), "none")) %>%
    arrange(counts_per_cov_plot)
  sorted_counts\$ids <- seq_len(nrow(sorted_counts))

  Q3 <- quantile(sorted_counts\$counts_per_cov_plot, 0.90, na.rm = TRUE)
  if (n_no_count / n_total > 0.10) {
    Q1 <- NA
    logdiff_text <- "LogDiff: not defined (more than 10 % of variants without counts)"
  } else {
    Q1 <- quantile(sorted_counts\$counts_per_cov_plot, 0.10, na.rm = TRUE)
    logdiff_text <- paste0("LogDiff = ", round(log10(Q3) - log10(Q1), 2))
  }
  percentile_lines <- c(Q1, Q3)[!is.na(c(Q1, Q3))]

  label <- paste0(
    logdiff_text, "\\n",
    "No counts: ", n_no_count, " of ", n_total, " variants (", round(100 * n_no_count / n_total, 1), " %)",
    if (n_zeroed > 0) paste0(", ", n_zeroed, " of them corrected to 0") else ""
  )
  caption <- paste0(
    "Variants without counts are drawn at half the lowest observed counts per coverage (",
    format(signif(floor_value, 3), scientific = TRUE), ", dashed line)."
  )
  y_top <- max(sorted_counts\$counts_per_cov_plot)

  line_plot <- ggplot(sorted_counts, aes(x = ids, y = counts_per_cov_plot)) +
    geom_line(color = "black") +
    geom_hline(yintercept = floor_value, linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = percentile_lines, linetype = "dotted", color = "black") +
    annotate("text", x = 0, y = y_top, label = label, hjust = 0, vjust = 1, size = 4.5, color = "black") +
    scale_y_continuous(trans = "log10") +
    labs(x = "Variants", y = "Counts per Coverage", caption = caption) +
    theme_minimal()
  ggsave(filename = "logdiff_plot.pdf", plot = line_plot, device = "pdf", width = 10, height = 6)

  colored_plot <- ggplot(sorted_counts, aes(x = ids, y = counts_per_cov_plot, color = category)) +
    geom_hline(yintercept = floor_value, linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = percentile_lines, linetype = "dotted", color = "black") +
    geom_point(size = 0.9) +
    annotate("text", x = 1, y = y_top, label = label, hjust = 0, vjust = 1, size = 4.5, color = "black") +
    scale_y_continuous(trans = "log10") +
    scale_color_manual(values = c("1" = "chocolate", "2" = "darkolivegreen3", "3" = "deepskyblue1", "none" = "grey70"),
                       breaks = c("1", "2", "3", "none"),
                       labels = c("1", "2", "3", "no counts"),
                       name = "Varying \\n Bases") +
    labs(x = "Variants", y = "Counts per Coverage", caption = caption) +
    theme_minimal()
  ggsave(filename = "logdiff_varying_bases.pdf", plot = colored_plot, device = "pdf", width = 10, height = 6)
}

#####
# run function
#####
logdiff_plot(
  completed_path = "$library_completed_variantCounts"
)

####
# create versions.yml
####
r_version <- strsplit(version[['version.string']], ' ')[[1]][3]
dplyr_version <- as.character(packageVersion("dplyr"))
ggplot2_version <- as.character(packageVersion("ggplot2"))
scales_version <- as.character(packageVersion("scales"))

if (is.null(r_version)) r_version <- "unknown"
if (length(dplyr_version) == 0) dplyr_version <- "unknown"
if (length(ggplot2_version) == 0) ggplot2_version <- "unknown"
if (length(scales_version) == 0) scales_version <- "unknown"

f <- file("versions.yml", "w")
writeLines(
  c(
    '"${task.process}":',
    paste('    r-base:', r_version),
    paste('    r-dplyr:', dplyr_version),
    paste('    r-ggplot2:', ggplot2_version),
    paste('    r-scales:', scales_version)
  ),
  f
)
close(f)
