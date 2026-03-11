#!/usr/bin/env Rscript

## fitness QC plots for nf-core/deepmutscan
## 28.10.2025
## maximilian.stammnitz@crg.eu

## lower panels: scatter + x=y (log-log version for counts)
panel_xy_abline_counts <- function(x, y, ...) {
  op <- par("xpd"); on.exit(par(xpd = op), add = TRUE)
  par(xpd = FALSE)
  points(x, y, pch = 16, cex = 0.1, ...)
  abline(a = 0, b = 1, lty = 2, col = "grey50")
}

## upper panels: Pearson r (log-log-transformed)
panel_cor_counts <- function(x, y, digits = 2, prefix = "r = ", cex.text = 1.4, ...) {
  r <- suppressWarnings(cor(log(x), log(y), use = "pairwise.complete.obs", method = "pearson"))
  lab <- if (is.finite(r)) bquote(italic(r) == .(round(r, digits))) else bquote(italic(r) == NA)

  ## Save/restore full graphics state we touch
  op <- par(c("usr", "xpd", "xlog", "ylog"))
  on.exit(par(op), add = TRUE)

  ## Draw in normalized 0..1 panel coords with logs OFF so text is visible
  par(xlog = FALSE, ylog = FALSE, xpd = FALSE, usr = c(0, 1, 0, 1))
  text(0.5, 0.5, labels = lab, cex = cex.text, font = 1, col = "black")
}

## lower panels: scatter + x=y (linear version for fitness)
panel_xy_abline_fitness <- function(x, y, ...) {
  op <- par("xpd"); on.exit(par(xpd = op), add = TRUE)
  par(xpd = FALSE)
  points(x, y, pch = 16, cex = 0.5, ...)
  abline(a = 0, b = 1, lty = 2, col = "grey50")
}

## upper panels: Pearson r (linear)
panel_cor_fitness <- function(x, y, digits = 2, prefix = "r = ", cex.text = 1.4, ...) {
  r <- suppressWarnings(cor(x, y, use = "pairwise.complete.obs", method = "pearson"))
  lab <- if (is.finite(r)) bquote(italic(r) == .(round(r, digits))) else bquote(italic(r) == NA)

  ## Save/restore full graphics state we touch
  op <- par(c("usr", "xpd", "xlog", "ylog"))
  on.exit(par(op), add = TRUE)

  ## Draw in normalized 0..1 panel coords with logs OFF so text is visible
  par(xlog = FALSE, ylog = FALSE, xpd = FALSE, usr = c(0, 1, 0, 1))
  text(0.5, 0.5, labels = lab, cex = cex.text, font = 1, col = "black")
}

#' Plot input/output count correlations and fitness replicate correlations
#'
#' @param fitness_table_path Path to the input table (fitness_estimation.tsv)
#' @param out_counts_corr_pdf Path to write the counts correlation PDF
#' @param out_fitness_corr_pdf Path to write the fitness correlation PDF
#'
#' @return Invisibly returns TRUE; writes the two PDFs.
run_fitness_plots <- function(fitness_table_path,
                              out_counts_corr_pdf,
                              out_fitness_corr_pdf) {

  merged.counts.fitness <- read.table(fitness_table_path, sep = "\\t", header = TRUE, check.names = FALSE)

  ## identify the right samples
  inputs  <- grep("input",  colnames(merged.counts.fitness))
  outputs <- grep("output", colnames(merged.counts.fitness))

  ## 5. Plot input vs. output counts ##
  #####################################
  pdf(out_counts_corr_pdf, height = 9, width = 14)
  pairs(merged.counts.fitness[, c(inputs, outputs)] + 1,  ## use a pseudo-count
        lower.panel = panel_xy_abline_counts,
        upper.panel = panel_cor_counts,
        cex.text = 2,
        log = "xy")
  dev.off()

  ## 6. Plot fitness correlations ##
  ##################################
  fitness.repl <- grep("rescaled_fitness", colnames(merged.counts.fitness))

  if (length(fitness.repl) > 1) {
    pdf(out_fitness_corr_pdf, height = 9, width = 14)
    pairs(merged.counts.fitness[, fitness.repl],
          lower.panel = panel_xy_abline_fitness,
          upper.panel = panel_cor_fitness,
          cex.text = 2,
          xlim = c(-3, 1),
          ylim = c(-3, 1))
    dev.off()
  } else {
    ## If only one (or zero) rescaled_fitness columns exist, still create an empty placeholder
    ## so Nextflow finds the declared output.
    pdf(out_fitness_corr_pdf, height = 9, width = 14)
    plot.new()
    title("No replicate fitness columns found (need ≥2 'rescaled_fitness...')\\nCreated placeholder PDF.")
    dev.off()
  }

  invisible(TRUE)
}


#####
# run function
#####
run_fitness_plots(
  fitness_table_path = "$fitness_estimation_tsv",
  out_counts_corr_pdf = "fitness_estimation_count_correlation.pdf",
  out_fitness_corr_pdf = "fitness_estimation_fitness_correlation.pdf"
)

####
# create versions.yml
####
r_version <- strsplit(version[['version.string']], ' ')[[1]][3]

if (is.null(r_version)) r_version <- "unknown"

f <- file("versions.yml", "w")
writeLines(
  c(
    '"${task.process}":',
    paste('    r-base:', r_version)
  ),
  f
)
close(f)
