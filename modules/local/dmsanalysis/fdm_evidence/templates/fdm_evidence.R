#!/usr/bin/env Rscript
## nf-core/deepmutscan - false-double-mutant (FDM) error evidence, one sequencing library.
##
## For every observed single-nucleotide library variant (SNV) j: F_j, the reads carrying j together with a
## programmed multi-nucleotide codon variant (an anchor) within max_fd_dist codons but not closer than min_fd_dist,
## and M_j, the anchor-bearing reads that span codon j (the exposure). Also the genuine double mutants
## (two programmed multi-nucleotide variants in one read) by codon distance, for the estimate step.
##
## The model and its code are the fdmcorrect R package (bioconda::dms-fdmcorrect,
## https://github.com/MaximilianStammnitz/dms-fdmcorrect); its documentation describes each step.
##
## Nextflow template: an R `\$` is a literal dollar; backslashes meant for R are doubled.

suppressPackageStartupMessages(library(fdmcorrect))

fdm_evidence(
    raw_counts       = "$raw_counts",
    library_variants = "$filtered",
    reference        = "$wt_fasta",
    id               = "${meta.id}",
    max_fd_dist      = ${max_fd_dist},
    min_fd_dist      = ${min_fd_dist},
    outdir           = "."
)

## --- versions ---
r_version <- strsplit(version[['version.string']], ' ')[[1]][3]
if (is.null(r_version)) r_version <- "unknown"
writeLines(c('"${task.process}":', paste("    r-base:", r_version),
             paste("    dms-fdmcorrect:", as.character(packageVersion("fdmcorrect")))), "versions.yml")
