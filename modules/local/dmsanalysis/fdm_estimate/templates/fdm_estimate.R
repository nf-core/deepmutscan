#!/usr/bin/env Rscript
## nf-core/deepmutscan - false-double-mutant (FDM) error rates, one pool of libraries.
##
## Sums the evidence of the pool's libraries per SNV and estimates every SNV's error rate e_j once:
## `eb`, a Gamma-Poisson empirical-Bayes estimate with hierarchical priors (global, substitution type and,
## with `trinucleotide_context`, the 96 contexts), or `mle`, F_j / M_j. With genuine-double subtraction,
## the reads expected from real double mutants are removed from F first. A pool needs at least 100
## false-double reads (`eb`) or 1,000 (`mle`); with fewer it is left uncorrected. A wildtype library is
## summarised without rates. The libraries are combined in the order of their names, whatever order
## their tasks finished in.
##
## The model and its code are the fdmcorrect R package (bioconda::dms-fdmcorrect,
## https://github.com/MaximilianStammnitz/dms-fdmcorrect); its documentation describes each step.
##
## Nextflow template: an R `\$` is a literal dollar; backslashes meant for R are doubled.

suppressPackageStartupMessages(library(fdmcorrect))

files <- function(x) Filter(nzchar, strsplit(x, " ", fixed = TRUE)[[1]])

fdm_estimate(
    evidence                 = files("${[evidence].flatten().join(' ')}"),
    profiles                 = files("${[profiles].flatten().join(' ')}"),
    exposures                = files("${[exposures].flatten().join(' ')}"),
    id                       = "${meta.id}",
    method                   = "${method}",
    prior                    = "${prior}",
    max_fd_dist              = ${max_fd_dist},
    subtract_genuine_doubles = "${subtract_genuine_doubles}" == "true",
    wildtype                 = "${meta.wildtype}" == "true",
    outdir                   = "."
)

## --- versions ---
r_version <- strsplit(version[['version.string']], ' ')[[1]][3]
if (is.null(r_version)) r_version <- "unknown"
writeLines(c('"${task.process}":', paste("    r-base:", r_version),
             paste("    dms-fdmcorrect:", as.character(packageVersion("fdmcorrect")))), "versions.yml")
