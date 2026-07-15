#!/usr/bin/env Rscript

## default fitness estimation for nf-core/deepmutscan
## 18.02.2026
## maximilian.stammnitz@crg.eu

## 0. Libraries ##
##################

suppressPackageStartupMessages({
  library(Biostrings)
})

## --- Helper functions ---

# calculate nt hamming distances from the specified WT
compute_nt_hamming <- function(merged.counts, wt.seq) {
  merged.counts <- cbind("nt_ham" = rep(NA, nrow(merged.counts)), merged.counts)
  for (i in 1:nrow(merged.counts)){
    tmp.wt <- strsplit(as.character(wt.seq), "")[[1]]
    tmp.mut <- strsplit(as.character(merged.counts\$nt_seq[i]), "")[[1]]
    if(length(which(tmp.mut != tmp.wt)) == 0){
      merged.counts\$nt_ham[i] <- 0
      rm(tmp.mut, tmp.wt)
      next
    }else{
      merged.counts\$nt_ham[i] <- length(which(tmp.mut != tmp.wt))
      rm(tmp.mut, tmp.wt)
      next
    }
  }
  merged.counts
}

# translate sequences and add aa_seq
add_aa_seq <- function(merged.counts) {
  merged.counts <- cbind("aa_seq" = as.character(translate(DNAStringSet(merged.counts\$nt_seq))), merged.counts)
  merged.counts
}

# calculate AA hamming distances from the WT
compute_aa_hamming <- function(merged.counts, wt.seq.aa) {
  merged.counts <- cbind("aa_ham" = rep(NA, nrow(merged.counts)), merged.counts)
  for (i in 1:nrow(merged.counts)){
    tmp.wt <- strsplit(as.character(wt.seq.aa), "")[[1]]
    tmp.mut <- strsplit(as.character(merged.counts\$aa_seq[i]), "")[[1]]
    if(length(which(tmp.mut != tmp.wt)) == 0){
      merged.counts\$aa_ham[i] <- 0
      rm(tmp.mut, tmp.wt)
      next
    }else{
      merged.counts\$aa_ham[i] <- length(which(tmp.mut != tmp.wt))
      rm(tmp.mut, tmp.wt)
      next
    }
  }
  merged.counts
}

# name the mutations
name_mutations <- function(merged.counts, wt.seq.aa) {
  merged.counts <- cbind("wt_aa" = rep(NA, nrow(merged.counts)),
                         "pos" = rep(NA, nrow(merged.counts)),
                         "mut_aa" = rep(NA, nrow(merged.counts)), merged.counts)
  for (i in 1:nrow(merged.counts)){
    if(merged.counts\$aa_ham[i] == 0){
      next
    }else{
      tmp.wt <- strsplit(as.character(wt.seq.aa), "")[[1]]
      tmp.mut <- strsplit(as.character(merged.counts\$aa_seq[i]), "")[[1]]
      merged.counts\$pos[i] <- which(tmp.mut != tmp.wt)
      merged.counts\$'wt_aa'[i] <- tmp.wt[merged.counts\$pos[i]]
      merged.counts\$'mut_aa'[i] <- tmp.mut[merged.counts\$pos[i]]
      rm(tmp.mut, tmp.wt)
    }
  }
  merged.counts
}

# find stops, WT and WT; aggregate counts of variants which are identical on the aa (but not nt) level
aggregate_by_aa <- function(merged.counts) {
  ## find stops, WT and WT
  merged.counts <- cbind(merged.counts,
                         "wt" = rep(NA, nrow(merged.counts)),
                         "stop" = rep(NA, nrow(merged.counts)))
  merged.counts\$wt[which(merged.counts\$nt_ham == 0)] <- TRUE
  merged.counts\$stop[which(merged.counts\$'mut_aa' == "*")] <- TRUE

  ## aggregate counts of variants which are identical on the aa (but not nt) level
  ## exception: wildtype ones
  ## thereby shrinking the matrix
  uniq.aa.vars <- unique(merged.counts\$aa_seq)
  uniq.aa.vars <- uniq.aa.vars[-which(uniq.aa.vars == merged.counts\$aa_seq[which(merged.counts\$wt == TRUE)])]
  for(i in 1:length(uniq.aa.vars)){
    tmp.aa_seq <- uniq.aa.vars[i]
    hits <- which(as.character(merged.counts\$aa_seq) == tmp.aa_seq)
    if(length(hits) == 1){
      rm(tmp.aa_seq, hits)
      next
    }else{
      for(j in grep("input|output", colnames(merged.counts))){
        merged.counts[hits[1],j] <- sum(merged.counts[hits,j], na.rm = TRUE)
      }
      merged.counts[hits[1], "nt_seq"] <- paste(merged.counts[hits, "nt_seq"], collapse = ", ")
      merged.counts <- merged.counts[-hits[-1],]
      rm(tmp.aa_seq, hits)
      next
    }
  }
  merged.counts
}

# calculate bimodal fitness distribution
density_peaks <- function(x, adjust = 1, ...) {

  ## obtain density distribution
  d <- density(x, adjust = adjust, n = 5000, na.rm = T, ...)
  y <- d\$y
  idx <- which(diff(sign(diff(y))) == -2) + 1
  if (length(idx) == 0) return(numeric(0))

  ## order density peaks by height (y value); take top 2
  idx <- idx[order(y[idx], decreasing = TRUE)]
  peaks_x <- d\$x[idx]
  peaks_x[seq_len(min(2, length(peaks_x)))]

}

# 3. Raw fitness calculations ##
calc_raw_fitness <- function(merged.counts, exp.design) {
  ## how many fitness replicates are there
  reps <- length(unique(exp.design\$experiment_replicate))
  for (i in 1:reps){
    merged.counts <- cbind(merged.counts, rep(NA, nrow(merged.counts)))
    colnames(merged.counts)[ncol(merged.counts)] <- paste0("raw_fitness_rep", i)
  }

  ## calculate raw fitness of all variants vs. WT variant
  for (i in 1:reps){

    ### collect counts
    tmp.input.counts <- merged.counts[,paste0("input", i)]
    tmp.output.counts <- merged.counts[,paste0("output", i)]

    ### add pseudo-count to zero-outputs (if the corresponding input count is non-zero)
    tmp.output.counts[which(tmp.output.counts == 0 & tmp.input.counts != 0)] <- 1

    ### take logs
    tmp.wt.log.ratio <- log(tmp.output.counts[which(merged.counts\$wt == TRUE)] /
                              tmp.input.counts[which(merged.counts\$wt == TRUE)])
    tmp.fitness <- log(tmp.output.counts /
                         tmp.input.counts) - tmp.wt.log.ratio

    ### uncertain values to NA
    tmp.fitness[which(is.na(tmp.fitness) == TRUE)] <- NA
    tmp.fitness[which(tmp.fitness == "Inf")] <- NA

    ### add to table
    merged.counts[,c(ncol(merged.counts) - reps + i)] <- tmp.fitness

    ### clean up
    rm(tmp.fitness, tmp.wt.log.ratio, tmp.output.counts, tmp.input.counts)
  }

  list(merged.counts = merged.counts, reps = reps)
}

# 4. Fitness and error refinements ##
rescale_and_summarize <- function(merged.counts, reps) {
  ## center the raw fitness distributions on 0 (median of wildtype synonymous) and -1 (median of stops)
  for (i in 1:reps){

    merged.counts <- cbind(merged.counts, rep(NA, nrow(merged.counts)))
    colnames(merged.counts)[ncol(merged.counts)] <- paste0("rescaled_fitness_rep", i)

    ### anchors: median fitness of synonymous (aa_ham == 0) -> 0, of stops -> -1
    raw_col <- ncol(merged.counts) - reps
    tmp.raw <- merged.counts[,raw_col]
    tmp.wt.fitness.med   <- median(tmp.raw[which(merged.counts\$aa_ham == 0)], na.rm = TRUE)
    tmp.stop.fitness.med <- median(tmp.raw[which(merged.counts\$stop == TRUE)], na.rm = TRUE)

    ### if an anchor is missing, or the two anchors coincide (e.g. on sparse data,
    ### where lm() would return an NA slope), fall back to the modes of the
    ### (bimodal) fitness distribution to define the missing/degenerate anchor
    if(is.na(tmp.wt.fitness.med) | is.na(tmp.stop.fitness.med) |
       (!is.na(tmp.wt.fitness.med) & !is.na(tmp.stop.fitness.med) & tmp.wt.fitness.med == tmp.stop.fitness.med)){
      tmp.peaks <- sort(density_peaks(x = tmp.raw))
      tmp.hi <- if(!is.na(tmp.wt.fitness.med)) tmp.wt.fitness.med else if(length(tmp.peaks) >= 1) tmp.peaks[length(tmp.peaks)] else NA
      if(!is.na(tmp.stop.fitness.med) & (is.na(tmp.hi) || tmp.stop.fitness.med != tmp.hi)){
        tmp.lo <- tmp.stop.fitness.med
      }else{
        tmp.cand <- tmp.peaks[!is.na(tmp.peaks) & tmp.peaks != tmp.hi]
        tmp.lo <- if(length(tmp.cand) >= 1) tmp.cand[1] else NA
      }
      tmp.wt.fitness.med <- tmp.hi
      tmp.stop.fitness.med <- tmp.lo
      rm(tmp.peaks, tmp.hi)
    }

    ### two distinct anchors -> linear rescale; a single usable anchor -> WT-centre only
    if(!is.na(tmp.wt.fitness.med) & !is.na(tmp.stop.fitness.med) & tmp.wt.fitness.med != tmp.stop.fitness.med){
      lm.rescale <- lm(c(0, -1) ~ c(tmp.wt.fitness.med, tmp.stop.fitness.med))
      merged.counts[,ncol(merged.counts)] <- tmp.raw * lm.rescale\$coefficients[[2]] + lm.rescale\$coefficients[[1]]
      rm(lm.rescale)
    }else if(!is.na(tmp.wt.fitness.med)){
      merged.counts[,ncol(merged.counts)] <- tmp.raw - tmp.wt.fitness.med
    }

    rm(raw_col, tmp.raw, tmp.wt.fitness.med, tmp.stop.fitness.med)
  }

  ## calculate fitness mean and standard deviation across replicates
  merged.counts <- cbind(merged.counts,
                         "mean fitness" = rep(NA, nrow(merged.counts)),
                         "fitness sd" = rep(NA, nrow(merged.counts)))

  if(reps == 1){

    merged.counts\$'mean fitness' <- merged.counts[,ncol(merged.counts) - 2]

  }else if(reps > 1){

    merged.counts\$'mean fitness' <- apply(merged.counts[,c(ncol(merged.counts) - 1 - reps):c(ncol(merged.counts) - 2)],
                                          1,
                                          mean,
                                          na.rm = TRUE)
    merged.counts\$'fitness sd' <- apply(merged.counts[,c(ncol(merged.counts) - 1 - reps):c(ncol(merged.counts) - 2)],
                                        1,
                                        sd,
                                        na.rm = TRUE)

  }

  merged.counts
}

## --- Main function ---

#' Run default fitness estimation with configurable I/O paths
#'
#' @param counts_path Path to counts_merged.tsv
#' @param design_path Path to experimentalDesign.tsv
#' @param wt_seq_path Path to synonymous_wt.txt (single line DNA sequence)
#' @param output_path Path to write fitness_estimation.tsv
#'
#' @return Invisibly returns the final data.frame; writes the output to output_path.
run_fitness_estimation <- function(counts_path,
                                   design_path,
                                   wt_seq_path,
                                   output_path) {
  ## 1. Import key files ##
  #########################

  merged.counts <- read.table(counts_path, sep = "\t", header = TRUE, check.names = FALSE)
  exp.design   <- read.table(design_path, sep = "\t", header = TRUE, check.names = FALSE)
  wt.seq       <- DNAString(as.character(read.table(wt_seq_path)))
  wt.seq.aa    <- translate(wt.seq)

  ## 2. Pre-processing the count table ##
  #######################################

  ## calculate nt hamming distances from the specified WT
  merged.counts <- compute_nt_hamming(merged.counts, wt.seq)

  ## translate sequences
  merged.counts <- add_aa_seq(merged.counts)

  ## calculate AA hamming distances from the WT
  merged.counts <- compute_aa_hamming(merged.counts, wt.seq.aa)

  ## name the mutations
  merged.counts <- name_mutations(merged.counts, wt.seq.aa)

  ## find stops, WT and WT; aggregate AA-identical variants (except WT)
  merged.counts <- aggregate_by_aa(merged.counts)

  ## 3. Raw fitness calculations ##
  #################################
  fitness_res <- calc_raw_fitness(merged.counts, exp.design)
  merged.counts <- fitness_res\$merged.counts
  reps <- fitness_res\$reps

  ## 4. Fitness and error refinements ##
  ######################################
  merged.counts <- rescale_and_summarize(merged.counts, reps)

  ## clean up
  rm(reps)

  ## export
  write.table(merged.counts, output_path,
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t", na = "")

  invisible(merged.counts)
}


## 5. Version ##
################

# sessionInfo()
# R version 4.5.1 (2025-06-13)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sonoma 14.6.1
#
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#
# time zone: Europe/Madrid
# tzcode source: internal
#
# attached base packages:
# [1] stats4    stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
# [1] Biostrings_2.76.0   GenomeInfoDb_1.44.2 XVector_0.48.0      IRanges_2.42.0      S4Vectors_0.46.0
# [6] BiocGenerics_0.54.0 generics_0.1.4
#
# loaded via a namespace (and not attached):
# [1] httr_1.4.7              compiler_4.5.1          R6_2.6.1                tools_4.5.1
# [5] GenomeInfoDbData_1.2.14 rstudioapi_0.17.1       crayon_1.5.3            UCSC.utils_1.4.0
# [9] jsonlite_2.0.0



#####
# run function
#####
run_fitness_estimation(
  counts_path = "$counts_merged",
  design_path = "$exp_design",
  wt_seq_path = "$syn_wt_txt",
  output_path = "fitness_estimation.tsv"
)

####
# create versions.yml
####
r_version <- strsplit(version[['version.string']], ' ')[[1]][3]
Biostrings_version <- as.character(packageVersion("Biostrings"))

if (is.null(r_version)) r_version <- "unknown"
if (length(Biostrings_version) == 0) Biostrings_version <- "unknown"

f <- file("versions.yml", "w")
writeLines(
  c(
    '"\${task.process}":',
    paste('    r-base:', r_version),
    paste('    r-Biostrings:', Biostrings_version)
  ),
  f
)
close(f)
