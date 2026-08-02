#!/usr/bin/env Rscript

## false-double mutant codon sequencing-error correction in nf-core/deepmutscan (nucleotide level).
##
## Two estimators, ported VERBATIM from Maxi's scripts (24.06.2026):
##   * MLE (default): per-SNV maximum-likelihood error rate e = F / m, applied per variant.
##   * EB:            pool F and m across SNVs, then a Gamma-Poisson empirical-Bayes fit with
##                    six-category shrinkage -> posterior-mean rate e_EB.
## `$method` (mle|eb) selects which one runs. Everything except I/O plumbing and the dispatch is
## kept exactly as written.
##
## Nextflow template: an R `\$` is a literal dollar; backslashes meant for R are doubled.

suppressMessages({
  library(Biostrings)
  library(stringr)
})

## -------------------------------------------------------------------------------------------------
## MLE estimator (verbatim)
## -------------------------------------------------------------------------------------------------
seq_error_correct_by_false_doubles_MLE <- function(wt_path, input_count_path_raw, input_count_path_processed, output_file_path, seq_error_rate_path, codon_window){

  ## load data (nucleotide-level counts from GATK)

  # WT sequence
  wt.seq <- readDNAStringSet(wt_path)

  # for the seq error rates, build a vector to be filled
  seq.error.names <- c()
  for(i in 1:nchar(wt.seq)){
    tmp.wt <- strsplit(as.character(wt.seq), "")[[1]][i]
    seq.error.names <- c(seq.error.names, paste0(i, ":", tmp.wt, ">", c("A", "T", "G", "C"))[-match(tmp.wt, c("A", "T", "G", "C"))])
  }
  seq.error.rate <- matrix(0, nrow = length(seq.error.names), ncol = 2)
  rownames(seq.error.rate) <- seq.error.names
  colnames(seq.error.rate) <- c("1nt counts/coverage", "1nt false counts/coverage")

  # Set the column names
  colnames <- c("counts", "cov", "mean_length_variant_reads", "varying_bases",
                "base_mut", "varying_codons", "codon_mut", "aa_mut", "pos_mut")
  input.counts.raw <- read.table(input_count_path_raw, sep = "\\t", header = FALSE, fill = TRUE, col.names = colnames)
  input.counts.processed <- read.csv(input_count_path_processed)

  ## append key columns
  input.counts.processed\$counts_corrected <- input.counts.processed\$counts
  input.counts.processed\$counts_per_cov_corrected <- input.counts.processed\$counts_per_cov

  ## calculate some things upfront

  ### highest accessible codon
  max.codon <- max(as.numeric(str_split_fixed(input.counts.raw\$codon_mut, ":", 2)[,1]), na.rm = T)

  ### total distribution of false double mutant coverage decay (%) vs. codon-codon distance
  cat("Estimating the distance-dependent coverage decay of false double mutants...\\n")
  all.false.doubles <- input.counts.raw
  all.false.doubles <- all.false.doubles[which(all.false.doubles[,"varying_bases"] >= 3 & all.false.doubles[,"varying_bases"] < 5 & all.false.doubles[,"varying_codons"] == 2),,drop = F]
  all.false.doubles\$codon_dist <- vapply(regmatches(all.false.doubles[,"codon_mut"], gregexpr("\\\\d+(?=:)", all.false.doubles[,"codon_mut"], perl = TRUE)), function(z) abs(diff(as.integer(z))), integer(1))
  
  ### additional filters: 

  #### remove outlier double codon mutants with very high counts
  all.false.doubles <- all.false.doubles[which(all.false.doubles\$counts <= quantile(all.false.doubles\$counts, c(0.995))),]

  #### only keep 2+1 nt and 3+1 nt double codon mutants
  all.false.doubles.nt <- all.false.doubles\$codon_mut
  all.false.doubles.nt <- str_split_fixed(all.false.doubles.nt, ", ", 2)
  all.false.doubles.nt[,1] <- str_split_fixed(all.false.doubles.nt[,1], ":", 2)[,2]
  all.false.doubles.nt[,2] <- str_split_fixed(all.false.doubles.nt[,2], ":", 2)[,2]
  classify_double_codon_variant <- function(x) {
    parts <- strsplit(x, ">", fixed = TRUE)
    n_changes <- vapply(parts, function(p) {
      if (length(p) != 2L || nchar(p[1]) != nchar(p[2])) {
        return(NA_integer_)
      }
      sum(strsplit(p[1], "")[[1]] != strsplit(p[2], "")[[1]])
    }, integer(1))
    labels <- c(1, 2, 3)
    ifelse(n_changes %in% 1:3,
           labels[n_changes],
           ifelse(n_changes == 0L, "no change", NA_character_))
  }
  all.false.doubles.nt[,1] <- classify_double_codon_variant(all.false.doubles.nt[,1])
  all.false.doubles.nt[,2] <- classify_double_codon_variant(all.false.doubles.nt[,2])
  if(length(which(all.false.doubles.nt[,1] == 2 & all.false.doubles.nt[,2] == 2) != 0)){
    all.false.doubles <- all.false.doubles[-which(all.false.doubles.nt[,1] == 2 & all.false.doubles.nt[,2] == 2),]
    all.false.doubles.nt <- all.false.doubles.nt[-which(all.false.doubles.nt[,1] == 2 & all.false.doubles.nt[,2] == 2),]
  }
  if(length(which(all.false.doubles.nt[,1] == 1 & all.false.doubles.nt[,2] == 1) != 0)){
    all.false.doubles <- all.false.doubles[-which(all.false.doubles.nt[,1] == 1 & all.false.doubles.nt[,2] == 1),]
    all.false.doubles.nt <- all.false.doubles.nt[-which(all.false.doubles.nt[,1] == 1 & all.false.doubles.nt[,2] == 1),]
  }
  if(length(which(c(all.false.doubles.nt[,1] == 3 & all.false.doubles.nt[,2] == 2) | c(all.false.doubles.nt[,1] == 2 & all.false.doubles.nt[,2] == 3)) != 0)){
    all.false.doubles <- all.false.doubles[-which(c(all.false.doubles.nt[,1] == 3 & all.false.doubles.nt[,2] == 2) | c(all.false.doubles.nt[,1] == 2 & all.false.doubles.nt[,2] == 3)),]
    all.false.doubles.nt <- all.false.doubles.nt[-which(c(all.false.doubles.nt[,1] == 3 & all.false.doubles.nt[,2] == 2) | c(all.false.doubles.nt[,1] == 2 & all.false.doubles.nt[,2] == 3)),]
  }
  if(length(which(is.na(all.false.doubles.nt[,1]) == T | is.na(all.false.doubles.nt[,2]) == T) != 0)){
    all.false.doubles <- all.false.doubles[-which(is.na(all.false.doubles.nt[,1]) == T | is.na(all.false.doubles.nt[,2]) == T),]
    all.false.doubles.nt <- all.false.doubles.nt[-which(is.na(all.false.doubles.nt[,1]) == T | is.na(all.false.doubles.nt[,2]) == T),]
  }

  all.false.doubles.codons <- str_split_fixed(all.false.doubles\$codon_mut, ", ", 2)
  all.false.doubles.codons[,1] <- as.integer(sub(":.*", "", all.false.doubles.codons[,1]))
  all.false.doubles.codons[,2] <- as.integer(sub(":.*", "", all.false.doubles.codons[,2]))
  median.high.conf.coverage.per.pos <- rep(NA, max(as.numeric(str_split_fixed(input.counts.processed\$codon_mut, ":", 2)[,1])))
  names(median.high.conf.coverage.per.pos) <- 1:max(as.numeric(str_split_fixed(input.counts.processed\$codon_mut, ":", 2)[,1]))
  for(i in 1:length(median.high.conf.coverage.per.pos)){
    median.high.conf.coverage.per.pos[i] <- median(input.counts.processed[which(names(median.high.conf.coverage.per.pos)[i] == str_split_fixed(input.counts.processed\$codon_mut, ":", 2)[,1]),"cov"])
  }
  all.false.doubles.codons[,1] <- median.high.conf.coverage.per.pos[match(all.false.doubles.codons[,1], names(median.high.conf.coverage.per.pos))]
  all.false.doubles.codons[,2] <- median.high.conf.coverage.per.pos[match(all.false.doubles.codons[,2], names(median.high.conf.coverage.per.pos))]
  all.false.doubles\$max_cov_possible <- as.numeric(apply(all.false.doubles.codons, 1, min))
  all.false.doubles\$perc_cov_to_max <- 100 * all.false.doubles\$cov / all.false.doubles\$max_cov_possible

  ## deepmutscan guard (not in the verbatim script): a run can legitimately contain NO observed false
  ## double mutants -- a small, shallow or tightly-filtered library, or the nf-core test dataset, where
  ## every variant is single-codon. The verbatim code then hands lm() an empty frame and dies with
  ## "0 (non-NA) cases". With no false doubles anywhere there is also no evidence of codon-level
  ## sequencing error to subtract, so the correct answer is a no-op: corrected == raw. The quadratic
  ## needs >= 3 distinct codon distances to be estimable at all, so require that too.
  .usable <- is.finite(log(all.false.doubles\$perc_cov_to_max))
  if (nrow(all.false.doubles) == 0 || length(unique(all.false.doubles\$codon_dist[.usable])) < 3) {
    cat("WARNING: only", sum(.usable), "usable false double mutant(s) found - the distance-dependent\\n",
        "coverage decay is not estimable, so NO sequencing-error correction is applied and the\\n",
        "corrected counts are identical to the raw counts.\\n")
    write.csv(seq.error.rate, file = seq_error_rate_path, row.names = T)
    write.csv(input.counts.processed, file = output_file_path, row.names = F)
    return(invisible(NULL))
  }

  fit <- lm(log(all.false.doubles\$perc_cov_to_max) ~ all.false.doubles\$codon_dist + I(all.false.doubles\$codon_dist^2))
  all.false.doubles\$pred_cov_perc <- exp(predict(fit))

  ### convert this into a look-up table: "% max. possible coverage" depending on codon-codon distance
  dual.codon.coverage.estimate <- matrix(NA, nrow = max(all.false.doubles\$codon_dist), ncol = 2)
  colnames(dual.codon.coverage.estimate) <- c("Codon distance", "Estimate % max. possible coverage")
  dual.codon.coverage.estimate[,"Codon distance"] <- 1:nrow(dual.codon.coverage.estimate)
  ## deepmutscan fix (Maxi-approved deviation from the verbatim script): fill the coverage-decay
  ## lookup from the fitted quadratic evaluated at EVERY integer codon distance, instead of match()ing
  ## only the distances that happen to have an observed double. `pred_cov_perc = exp(predict(fit))`
  ## depends solely on codon_dist, so this is bit-identical to the original at every observed distance
  ## and merely fills the previously-NA gaps -- those NAs otherwise crash MLE and leave NA corrected
  ## counts in EB (a distance with no observed double, e.g. around codon 288 on the GID1A min_counts=10
  ## data).
  .dccd <- dual.codon.coverage.estimate[,"Codon distance"]
  .cf <- coef(fit)
  dual.codon.coverage.estimate[,"Estimate % max. possible coverage"] <- exp(.cf[1] + .cf[2] * .dccd + .cf[3] * .dccd^2)

  ## Process the GATK file, only look at single nucleotide variants
  cat("Sequencing error correction of GATK counts...\\n")
  for (i in grep("[,]", input.counts.processed[,"base_mut"], invert = T)){

    ### algorithm:
    ### 1.) isolate 1nt variant of interest
    ### 2.) look for 2/3nt hits in a pre-specified window of codons upstream or downstream
    ### -> if there are zero false double mutants: skip / set correction factor to zero
    ### -> if there are any false double mutants: continue
    ### 3.) generate a table for ALL possible 2/3nt variants in the selected window:
    ### true high-conf. variant count | true high-conf. variant coverage | false double variant count (often 0 or 1) | false double variant coverage
    ### 4.) Maximum likelihood estimation (MLE) over all events
    ### 5.) Apply correction factor

    ## 1.) isolate 1nt variant of interest
    tmp.single <- input.counts.processed[i,"base_mut"]
    seq.error.rate[tmp.single, "1nt counts/coverage"] <- input.counts.processed[i,"counts_per_cov"]

    ## 2.) look for 2/3nt hits in up to 40 codons (120 bp) upstream or downstream
    tmp.false.doubles <- all.false.doubles[grep(paste0("(?<!\\\\d)",input.counts.processed[i,"codon_mut"]), all.false.doubles[,"codon_mut"], perl = T),,drop = F]

    ## set and enforce distance threshold
    tmp.ref.codon <- as.numeric(strsplit(input.counts.processed[i,"codon_mut"] , ":")[[1]][1])
    tmp.min.from.ref <- max(c(1, tmp.ref.codon - c(codon_window - 1)))
    tmp.max.from.ref <- min(c(max.codon, tmp.ref.codon + c(codon_window - 1)))
    tmp.false.doubles <- tmp.false.doubles[which(tmp.false.doubles\$codon_dist < codon_window),,drop = F]

    ## 3.) generate a table for ALL possible 2/3nt variants in the selected window
    tmp.summary.table.entries <- input.counts.processed
    tmp.summary.table.entries <- tmp.summary.table.entries[which(as.numeric(str_split_fixed(tmp.summary.table.entries\$codon_mut, ":", 2)[,1]) != tmp.ref.codon),]
    tmp.summary.table.entries <- tmp.summary.table.entries[which(as.numeric(str_split_fixed(tmp.summary.table.entries\$codon_mut, ":", 2)[,1]) >= tmp.min.from.ref &
                                                                   as.numeric(str_split_fixed(tmp.summary.table.entries\$codon_mut, ":", 2)[,1]) <= tmp.max.from.ref),]
    tmp.summary.table.entries <- tmp.summary.table.entries[which(tmp.summary.table.entries\$varying_bases > 1),]

    ### fill the table
    tmp.summary.table <- matrix(NA, ncol = 5, nrow = nrow(tmp.summary.table.entries))
    colnames(tmp.summary.table) <- c("codon distance",
                                     "true high-conf. variant count",
                                     "true high-conf. variant coverage",
                                     "false double variant count",
                                     "false double variant coverage")
    rownames(tmp.summary.table) <- paste0(c(tmp.summary.table.entries\$codon_mut),", ", c(tmp.summary.table.entries\$base_mut))

    tmp.summary.table[,"codon distance"] <- abs(as.numeric(str_split_fixed(tmp.summary.table.entries\$codon_mut, ":", 2)[,1]) - tmp.ref.codon)
    tmp.summary.table[,"true high-conf. variant count"] <- tmp.summary.table.entries\$counts
    tmp.summary.table[,"true high-conf. variant coverage"] <- tmp.summary.table.entries\$cov
    tmp.summary.table[,"false double variant count"] <- 0

    ## input the actual hits
    tmp.double.codon <- str_split_fixed(tmp.false.doubles\$codon_mut, ", ", 2)
    tmp.match1 <- match(tmp.double.codon[,1], str_split_fixed(rownames(tmp.summary.table), ", ", 2)[,1])
    tmp.match1[is.na(tmp.match1)] <- 0
    tmp.match2 <- match(tmp.double.codon[,2], str_split_fixed(rownames(tmp.summary.table), ", ", 2)[,1])
    tmp.match2[is.na(tmp.match2)] <- 0
    tmp.match <- c(tmp.match1 + tmp.match2)

    ### remove those false doubles for which there are no matched true high confidence calls
    if(any(tmp.match == 0)){
      tmp.false.doubles <- tmp.false.doubles[-which(tmp.match == 0),,drop = F]
      tmp.match <- tmp.match[-which(tmp.match == 0)]
      if(nrow(tmp.false.doubles) == 0){
        next
      }
    }

    tmp.summary.table[tmp.match,"false double variant count"] <- tmp.false.doubles[,"counts"]
    tmp.summary.table[tmp.match,"false double variant coverage"] <- tmp.false.doubles[,"cov"]
    tmp.summary.table <- as.data.frame(tmp.summary.table)

    ## 4.) Maximum likelihood estimation (MLE) over all events, iterating over codon positions

    ### summarise data by position
    tmp.summary.table.pos <- matrix(NA, nrow = length(unique(str_split_fixed(rownames(tmp.summary.table), ">", 2)[,1])), ncol = 5)
    colnames(tmp.summary.table.pos) <- colnames(tmp.summary.table)
    rownames(tmp.summary.table.pos) <- unique(str_split_fixed(rownames(tmp.summary.table), ">", 2)[,1])
    tmp.summary.table.pos[,"codon distance"] <- abs(as.numeric(str_split_fixed(rownames(tmp.summary.table.pos), ":", 2)[,1]) - as.numeric(str_split_fixed(input.counts.processed[i,"codon_mut"], ":", 2)[,1]))
    for(j in 1:nrow(tmp.summary.table.pos)){
      tmp.name <- rownames(tmp.summary.table.pos)[j]
      tmp.summary.table.pos[j,"true high-conf. variant count"] <- sum(tmp.summary.table[grep(paste0("^",tmp.name), rownames(tmp.summary.table)),"true high-conf. variant count"])
      tmp.summary.table.pos[j,"true high-conf. variant coverage"] <- max(tmp.summary.table[grep(paste0("^",tmp.name), rownames(tmp.summary.table)),"true high-conf. variant coverage"])
      tmp.summary.table.pos[j,"false double variant count"] <- sum(tmp.summary.table[grep(paste0("^",tmp.name), rownames(tmp.summary.table)),"false double variant count"])
      tmp.summary.table.pos[j,"false double variant coverage"] <- max(tmp.summary.table[grep(paste0("^",tmp.name), rownames(tmp.summary.table)),"false double variant coverage"], na.rm = T)

      ## infer a good proxy for the false double variant coverage with 0 counts (based estimated distance-dependent % coverage decay observed across all observed double mutants)
      if(tmp.summary.table.pos[j,"false double variant coverage"] == "-Inf"){
        tmp.summary.table.pos[j,"false double variant coverage"] <- c(dual.codon.coverage.estimate[abs(tmp.summary.table.pos[j,"codon distance"]),"Estimate % max. possible coverage"] / 100) * tmp.summary.table.pos[j,"true high-conf. variant coverage"]
        tmp.summary.table.pos[j,"false double variant coverage"] <- round(tmp.summary.table.pos[j,"false double variant coverage"])
      }
    }
    tmp.summary.table.pos <- as.data.frame(tmp.summary.table.pos)
    e_MLE <- sum(tmp.summary.table.pos\$`false double variant count`) /
      sum(tmp.summary.table.pos\$`true high-conf. variant count` * c(tmp.summary.table.pos\$`false double variant coverage` / tmp.summary.table.pos\$`true high-conf. variant coverage`))
    seq.error.rate[tmp.single,"1nt false counts/coverage"] <- e_MLE
    input.counts.processed[i,"counts_per_cov_corrected"] <- input.counts.processed[i,"counts_per_cov"] - e_MLE

    ## 5.) Apply correction factor
    input.counts.processed[i,"counts_corrected"] <- input.counts.processed[i,"counts"] - c(e_MLE * input.counts.processed[i,"cov"])

    ## round to nearest integer, do not allow for negative counts
    input.counts.processed[i,"counts_corrected"] <- round(input.counts.processed[i,"counts_corrected"])
    if(input.counts.processed[i,"counts_corrected"] < 0){
      input.counts.processed[i,"counts_corrected"] <- 0
      input.counts.processed[i,"counts_per_cov_corrected"] <- 0
    }

  }

  # Write the processed data
  write.csv(seq.error.rate, file = seq_error_rate_path, row.names = T)
  write.csv(input.counts.processed, file = output_file_path, row.names = F)

}

## -------------------------------------------------------------------------------------------------
## EB estimator (verbatim)
## -------------------------------------------------------------------------------------------------
seq_error_correct_by_false_doubles_EB <- function(wt_path, input_count_path_raw, input_count_path_processed, output_file_path, seq_error_rate_path, codon_window){

  ## load data (nucleotide-level counts from GATK)

  # WT sequence
  wt.seq <- readDNAStringSet(wt_path)

  # for the seq error rates, build a vector to be filled
  seq.error.names <- c()
  for(i in 1:nchar(wt.seq)){
    tmp.wt <- strsplit(as.character(wt.seq), "")[[1]][i]
    seq.error.names <- c(seq.error.names, paste0(i, ":", tmp.wt, ">", c("A", "T", "G", "C"))[-match(tmp.wt, c("A", "T", "G", "C"))])
  }
  seq.error.rate <- matrix(0, nrow = length(seq.error.names), ncol = 2)
  rownames(seq.error.rate) <- seq.error.names
  colnames(seq.error.rate) <- c("1nt counts/coverage", "1nt false counts/coverage")

  # Set the column names
  colnames <- c("counts", "cov", "mean_length_variant_reads", "varying_bases",
                "base_mut", "varying_codons", "codon_mut", "aa_mut", "pos_mut")
  input.counts.raw <- read.table(input_count_path_raw, sep = "\\t", header = FALSE, fill = TRUE, col.names = colnames)
  input.counts.processed <- read.csv(input_count_path_processed)

  ## append key columns
  input.counts.processed\$counts_corrected <- input.counts.processed\$counts
  input.counts.processed\$counts_per_cov_corrected <- input.counts.processed\$counts_per_cov

  ## calculate some things upfront

  ### highest accessible codon
  max.codon <- max(as.numeric(str_split_fixed(input.counts.raw\$codon_mut, ":", 2)[,1]), na.rm = T)

  ### total distribution of false double mutant coverage decay (%) vs. codon-codon distance
  cat("Estimating the distance-dependent coverage decay of false double mutants...\\n")
  all.false.doubles <- input.counts.raw
  all.false.doubles <- all.false.doubles[which(all.false.doubles[,"varying_bases"] >= 3 & all.false.doubles[,"varying_bases"] < 5 & all.false.doubles[,"varying_codons"] == 2),,drop = F]
  all.false.doubles\$codon_dist <- vapply(regmatches(all.false.doubles[,"codon_mut"], gregexpr("\\\\d+(?=:)", all.false.doubles[,"codon_mut"], perl = TRUE)), function(z) abs(diff(as.integer(z))), integer(1))
  
  ### additional filters: 

  #### remove outlier double codon mutants with very high counts
  all.false.doubles <- all.false.doubles[which(all.false.doubles\$counts <= quantile(all.false.doubles\$counts, c(0.995))),]

  #### only keep 2+1 nt and 3+1 nt double codon mutants
  all.false.doubles.nt <- all.false.doubles\$codon_mut
  all.false.doubles.nt <- str_split_fixed(all.false.doubles.nt, ", ", 2)
  all.false.doubles.nt[,1] <- str_split_fixed(all.false.doubles.nt[,1], ":", 2)[,2]
  all.false.doubles.nt[,2] <- str_split_fixed(all.false.doubles.nt[,2], ":", 2)[,2]
  classify_double_codon_variant <- function(x) {
    parts <- strsplit(x, ">", fixed = TRUE)
    n_changes <- vapply(parts, function(p) {
      if (length(p) != 2L || nchar(p[1]) != nchar(p[2])) {
        return(NA_integer_)
      }
      sum(strsplit(p[1], "")[[1]] != strsplit(p[2], "")[[1]])
    }, integer(1))
    labels <- c(1, 2, 3)
    ifelse(n_changes %in% 1:3,
           labels[n_changes],
           ifelse(n_changes == 0L, "no change", NA_character_))
  }
  all.false.doubles.nt[,1] <- classify_double_codon_variant(all.false.doubles.nt[,1])
  all.false.doubles.nt[,2] <- classify_double_codon_variant(all.false.doubles.nt[,2])
  if(length(which(all.false.doubles.nt[,1] == 2 & all.false.doubles.nt[,2] == 2) != 0)){
    all.false.doubles <- all.false.doubles[-which(all.false.doubles.nt[,1] == 2 & all.false.doubles.nt[,2] == 2),]
    all.false.doubles.nt <- all.false.doubles.nt[-which(all.false.doubles.nt[,1] == 2 & all.false.doubles.nt[,2] == 2),]
  }
  if(length(which(all.false.doubles.nt[,1] == 1 & all.false.doubles.nt[,2] == 1) != 0)){
    all.false.doubles <- all.false.doubles[-which(all.false.doubles.nt[,1] == 1 & all.false.doubles.nt[,2] == 1),]
    all.false.doubles.nt <- all.false.doubles.nt[-which(all.false.doubles.nt[,1] == 1 & all.false.doubles.nt[,2] == 1),]
  }
  if(length(which(c(all.false.doubles.nt[,1] == 3 & all.false.doubles.nt[,2] == 2) | c(all.false.doubles.nt[,1] == 2 & all.false.doubles.nt[,2] == 3)) != 0)){
    all.false.doubles <- all.false.doubles[-which(c(all.false.doubles.nt[,1] == 3 & all.false.doubles.nt[,2] == 2) | c(all.false.doubles.nt[,1] == 2 & all.false.doubles.nt[,2] == 3)),]
    all.false.doubles.nt <- all.false.doubles.nt[-which(c(all.false.doubles.nt[,1] == 3 & all.false.doubles.nt[,2] == 2) | c(all.false.doubles.nt[,1] == 2 & all.false.doubles.nt[,2] == 3)),]
  }
  if(length(which(is.na(all.false.doubles.nt[,1]) == T | is.na(all.false.doubles.nt[,2]) == T) != 0)){
    all.false.doubles <- all.false.doubles[-which(is.na(all.false.doubles.nt[,1]) == T | is.na(all.false.doubles.nt[,2]) == T),]
    all.false.doubles.nt <- all.false.doubles.nt[-which(is.na(all.false.doubles.nt[,1]) == T | is.na(all.false.doubles.nt[,2]) == T),]
  }

  all.false.doubles.codons <- str_split_fixed(all.false.doubles\$codon_mut, ", ", 2)
  all.false.doubles.codons[,1] <- as.integer(sub(":.*", "", all.false.doubles.codons[,1]))
  all.false.doubles.codons[,2] <- as.integer(sub(":.*", "", all.false.doubles.codons[,2]))
  median.high.conf.coverage.per.pos <- rep(NA, max(as.numeric(str_split_fixed(input.counts.processed\$codon_mut, ":", 2)[,1])))
  names(median.high.conf.coverage.per.pos) <- 1:max(as.numeric(str_split_fixed(input.counts.processed\$codon_mut, ":", 2)[,1]))
  for(i in 1:length(median.high.conf.coverage.per.pos)){
    median.high.conf.coverage.per.pos[i] <- median(input.counts.processed[which(names(median.high.conf.coverage.per.pos)[i] == str_split_fixed(input.counts.processed\$codon_mut, ":", 2)[,1]),"cov"])
  }
  all.false.doubles.codons[,1] <- median.high.conf.coverage.per.pos[match(all.false.doubles.codons[,1], names(median.high.conf.coverage.per.pos))]
  all.false.doubles.codons[,2] <- median.high.conf.coverage.per.pos[match(all.false.doubles.codons[,2], names(median.high.conf.coverage.per.pos))]
  all.false.doubles\$max_cov_possible <- as.numeric(apply(all.false.doubles.codons, 1, min))
  all.false.doubles\$perc_cov_to_max <- 100 * all.false.doubles\$cov / all.false.doubles\$max_cov_possible

  ## deepmutscan guard (not in the verbatim script) - see the identical guard in the MLE estimator.
  ## With no estimable coverage decay there is no evidence of codon-level sequencing error, so the
  ## correct answer is a no-op: corrected == raw, rather than lm() dying on an empty frame.
  .usable <- is.finite(log(all.false.doubles\$perc_cov_to_max))
  if (nrow(all.false.doubles) == 0 || length(unique(all.false.doubles\$codon_dist[.usable])) < 3) {
    cat("WARNING: only", sum(.usable), "usable false double mutant(s) found - the distance-dependent\\n",
        "coverage decay is not estimable, so NO sequencing-error correction is applied and the\\n",
        "corrected counts are identical to the raw counts.\\n")
    write.csv(seq.error.rate, file = seq_error_rate_path, row.names = T)
    write.csv(input.counts.processed, file = output_file_path, row.names = F)
    return(invisible(NULL))
  }

  fit <- lm(log(all.false.doubles\$perc_cov_to_max) ~ all.false.doubles\$codon_dist + I(all.false.doubles\$codon_dist^2))
  all.false.doubles\$pred_cov_perc <- exp(predict(fit))

  ### convert this into a look-up table: "% max. possible coverage" depending on codon-codon distance
  dual.codon.coverage.estimate <- matrix(NA, nrow = max(all.false.doubles\$codon_dist), ncol = 2)
  colnames(dual.codon.coverage.estimate) <- c("Codon distance", "Estimate % max. possible coverage")
  dual.codon.coverage.estimate[,"Codon distance"] <- 1:nrow(dual.codon.coverage.estimate)
  ## deepmutscan fix (Maxi-approved deviation from the verbatim script): fill the coverage-decay
  ## lookup from the fitted quadratic evaluated at EVERY integer codon distance, instead of match()ing
  ## only the distances that happen to have an observed double. `pred_cov_perc = exp(predict(fit))`
  ## depends solely on codon_dist, so this is bit-identical to the original at every observed distance
  ## and merely fills the previously-NA gaps -- those NAs otherwise crash MLE and leave NA corrected
  ## counts in EB (a distance with no observed double, e.g. around codon 288 on the GID1A min_counts=10
  ## data).
  .dccd <- dual.codon.coverage.estimate[,"Codon distance"]
  .cf <- coef(fit)
  dual.codon.coverage.estimate[,"Estimate % max. possible coverage"] <- exp(.cf[1] + .cf[2] * .dccd + .cf[3] * .dccd^2)

  ## Process the GATK file, only look at single nucleotide variants
  cat("Sequencing error correction of GATK counts...\\n")
  out.summary <- matrix(NA, ncol = 7, nrow = length(grep("[,]", input.counts.processed[,"base_mut"], invert = T)))
  colnames(out.summary) <- c("SNV", "Category",
                             "False_double_variant_count", "False_double_variant_coverage",
                             "True_high_conf_variant_count", "True_high_conf_variant_coverage",
                             "Exposure_m")
  out.summary[,"SNV"] <- input.counts.processed[grep("[,]", input.counts.processed[,"base_mut"], invert = T),"base_mut"]
  out.summary[,"Category"] <- str_split_fixed(out.summary[,"SNV"], ":", 2)[,2]
  category.translate <- cbind(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                              c("G>T", "G>C", "G>A", "A>T", "A>G", "A>C"))
  for(i in 1:nrow(out.summary)){
    if(out.summary[i,"Category"] %in% category.translate[,2]){
      out.summary[i,"Category"] <- category.translate[match(out.summary[i,"Category"],category.translate[,2]),1]
    }
  }

  for (i in grep("[,]", input.counts.processed[,"base_mut"], invert = T)){

    ### algorithm:
    ### 1.) isolate 1nt variant of interest
    ### 2.) look for 2/3nt hits in a pre-specified window of codons upstream or downstream
    ### -> if there are zero false double mutants: skip / set correction factor to zero
    ### -> if there are any false double mutants: continue
    ### 3.) generate a table for ALL possible 2/3nt variants in the selected window:
    ### true high-conf. variant count | true high-conf. variant coverage | false double variant count (often 0 or 1) | false double variant coverage
    ### 4.) Iterate of all positions
    ### 5.) Calculate the empirical Bayes
    ### 6.) Apply correction factor

    ## 1.) isolate 1nt variant of interest
    tmp.single <- input.counts.processed[i,"base_mut"]
    seq.error.rate[tmp.single, "1nt counts/coverage"] <- input.counts.processed[i,"counts_per_cov"]

    ## 2.) look for 2/3nt hits in up to 40 codons (120 bp) upstream or downstream
    tmp.false.doubles <- all.false.doubles[grep(paste0("(?<!\\\\d)",input.counts.processed[i,"codon_mut"]), all.false.doubles[,"codon_mut"], perl = T),,drop = F]
    tmp.ref.codon <- as.numeric(strsplit(input.counts.processed[i,"codon_mut"] , ":")[[1]][1])
    tmp.min.from.ref <- max(c(1, tmp.ref.codon - c(codon_window - 1)))
    tmp.max.from.ref <- min(c(max.codon, tmp.ref.codon + c(codon_window - 1)))
    tmp.false.doubles <- tmp.false.doubles[which(tmp.false.doubles\$codon_dist < codon_window),,drop = F]

    ## 3.) generate a table for ALL possible 2/3nt variants in the selected window
    tmp.summary.table.entries <- input.counts.processed
    tmp.summary.table.entries <- tmp.summary.table.entries[which(as.numeric(str_split_fixed(tmp.summary.table.entries\$codon_mut, ":", 2)[,1]) != tmp.ref.codon),]
    tmp.summary.table.entries <- tmp.summary.table.entries[which(as.numeric(str_split_fixed(tmp.summary.table.entries\$codon_mut, ":", 2)[,1]) >= tmp.min.from.ref &
                                                                   as.numeric(str_split_fixed(tmp.summary.table.entries\$codon_mut, ":", 2)[,1]) <= tmp.max.from.ref),]
    tmp.summary.table.entries <- tmp.summary.table.entries[which(tmp.summary.table.entries\$varying_bases > 1),]

    ### fill the table
    tmp.summary.table <- matrix(NA, ncol = 5, nrow = nrow(tmp.summary.table.entries))
    colnames(tmp.summary.table) <- c("codon distance",
                                     "true high-conf. variant count",
                                     "true high-conf. variant coverage",
                                     "false double variant count",
                                     "false double variant coverage")
    rownames(tmp.summary.table) <- paste0(c(tmp.summary.table.entries\$codon_mut),", ", c(tmp.summary.table.entries\$base_mut))

    tmp.summary.table[,"codon distance"] <- abs(as.numeric(str_split_fixed(tmp.summary.table.entries\$codon_mut, ":", 2)[,1]) - tmp.ref.codon)
    tmp.summary.table[,"true high-conf. variant count"] <- tmp.summary.table.entries\$counts
    tmp.summary.table[,"true high-conf. variant coverage"] <- tmp.summary.table.entries\$cov
    tmp.summary.table[,"false double variant count"] <- 0

    ## input the actual hits
    tmp.double.codon <- str_split_fixed(tmp.false.doubles\$codon_mut, ", ", 2)
    tmp.match1 <- match(tmp.double.codon[,1], str_split_fixed(rownames(tmp.summary.table), ", ", 2)[,1])
    tmp.match1[is.na(tmp.match1)] <- 0
    tmp.match2 <- match(tmp.double.codon[,2], str_split_fixed(rownames(tmp.summary.table), ", ", 2)[,1])
    tmp.match2[is.na(tmp.match2)] <- 0
    tmp.match <- c(tmp.match1 + tmp.match2)

    ### remove those false doubles for which there are no matched true high confidence calls
    if(any(tmp.match == 0)){
      tmp.false.doubles <- tmp.false.doubles[-which(tmp.match == 0),,drop = F]
      tmp.match <- tmp.match[-which(tmp.match == 0)]
      if(nrow(tmp.false.doubles) == 0){
        next
      }
    }

    tmp.summary.table[tmp.match,"false double variant count"] <- tmp.false.doubles[,"counts"]
    tmp.summary.table[tmp.match,"false double variant coverage"] <- tmp.false.doubles[,"cov"]
    tmp.summary.table <- as.data.frame(tmp.summary.table)

    ## 4.) Iterate over codon positions

    ### summarise data by position
    tmp.summary.table.pos <- matrix(NA, nrow = length(unique(str_split_fixed(rownames(tmp.summary.table), ">", 2)[,1])), ncol = 5)
    colnames(tmp.summary.table.pos) <- colnames(tmp.summary.table)
    rownames(tmp.summary.table.pos) <- unique(str_split_fixed(rownames(tmp.summary.table), ">", 2)[,1])
    tmp.summary.table.pos[,"codon distance"] <- abs(as.numeric(str_split_fixed(rownames(tmp.summary.table.pos), ":", 2)[,1]) - as.numeric(str_split_fixed(input.counts.processed[i,"codon_mut"], ":", 2)[,1]))
    for(j in 1:nrow(tmp.summary.table.pos)){
      tmp.name <- rownames(tmp.summary.table.pos)[j]
      tmp.summary.table.pos[j,"true high-conf. variant count"] <- sum(tmp.summary.table[grep(paste0("^",tmp.name), rownames(tmp.summary.table)),"true high-conf. variant count"])
      tmp.summary.table.pos[j,"true high-conf. variant coverage"] <- max(tmp.summary.table[grep(paste0("^",tmp.name), rownames(tmp.summary.table)),"true high-conf. variant coverage"])
      tmp.summary.table.pos[j,"false double variant count"] <- sum(tmp.summary.table[grep(paste0("^",tmp.name), rownames(tmp.summary.table)),"false double variant count"])
      tmp.summary.table.pos[j,"false double variant coverage"] <- max(tmp.summary.table[grep(paste0("^",tmp.name), rownames(tmp.summary.table)),"false double variant coverage"], na.rm = T)

      ## infer a good proxy for the false double variant coverage with 0 counts (based estimated distance-dependent % coverage decay observed across all observed double mutants)
      if(tmp.summary.table.pos[j,"false double variant coverage"] == "-Inf"){
        tmp.summary.table.pos[j,"false double variant coverage"] <- c(dual.codon.coverage.estimate[abs(tmp.summary.table.pos[j,"codon distance"]),"Estimate % max. possible coverage"] / 100) * tmp.summary.table.pos[j,"true high-conf. variant coverage"]
        tmp.summary.table.pos[j,"false double variant coverage"] <- round(tmp.summary.table.pos[j,"false double variant coverage"])
      }
    }
    tmp.summary.table.pos <- as.data.frame(tmp.summary.table.pos)

    ### fill the master table
    out.summary[match(tmp.single,out.summary[,"SNV"]),"False_double_variant_count"] <- sum(tmp.summary.table.pos\$`false double variant count`)
    out.summary[match(tmp.single,out.summary[,"SNV"]),"False_double_variant_coverage"] <- sum(tmp.summary.table.pos\$`false double variant coverage`)
    out.summary[match(tmp.single,out.summary[,"SNV"]),"True_high_conf_variant_count"] <- sum(tmp.summary.table.pos\$`true high-conf. variant count`)
    out.summary[match(tmp.single,out.summary[,"SNV"]),"True_high_conf_variant_coverage"] <- sum(tmp.summary.table.pos\$`true high-conf. variant coverage`)
    out.summary[match(tmp.single,out.summary[,"SNV"]),"Exposure_m"] <- sum(tmp.summary.table.pos\$`true high-conf. variant count` * c(tmp.summary.table.pos\$`false double variant coverage` / tmp.summary.table.pos\$`true high-conf. variant coverage`))

  }

  ## 5.) Calculate the empirical-Bayes sequencing error estimates
  out.summary <- out.summary[-which(is.na(out.summary[,"False_double_variant_count"]) == T),]
  out.summary <- as.data.frame(out.summary)
  class(out.summary\$False_double_variant_count) <- "numeric"
  class(out.summary\$False_double_variant_coverage) <- "numeric"
  class(out.summary\$True_high_conf_variant_count) <- "numeric"
  class(out.summary\$True_high_conf_variant_coverage) <- "numeric"
  class(out.summary\$Exposure_m) <- "numeric"

  ### calculate e_MLE across all
  out.summary\$e_MLE <- out.summary\$False_double_variant_count / out.summary\$Exposure_m

  ### calculate all the necessary rates, priors, etc.

  #### global background rate
  mu_global <- sum(out.summary\$False_double_variant_count) / sum(out.summary\$Exposure_m)
  mu_global <- max(mu_global, 1e-12) ##

  ##### marginal log-likelihood under Gamma-Poisson with fixed mean mu (ChatGPT function, adapted)
  ## prior: e ~ Gamma(shape = nu * mu, rate = nu)
  ## mean = mu, variance = mu / nu
  gp_loglik_nu <- function(log_nu, y, m, mu) {
    nu <- exp(log_nu)
    a <- max(nu * mu, 1e-12)

    # Negative-binomial marginal:
    # y ~ NB(size = a, prob = nu / (nu + m))
    ll <- lgamma(y + a) - lgamma(a) - lgamma(y + 1) +
      a * log(nu / (nu + m)) +
      y * log(m / (nu + m))

    -sum(ll)
  }
  opt_global <- optimize(f = gp_loglik_nu,
                         interval = c(log(1e-6), log(1e6)),
                         y = out.summary\$False_double_variant_count,
                         m = out.summary\$Exposure_m,
                         mu = mu_global)
  nu_global <- exp(opt_global\$minimum)

  #### category-level rates, shrunk to global rate
  mu_cat <- vector(mode = "list", length = 6)
  names(mu_cat) <- category.translate[,1]
  mu_cat_raw <- mu_cat
  for(i in 1:length(mu_cat)){
    mu_cat[[i]] <- out.summary[which(out.summary\$Category == names(mu_cat)[i]),]
    mu_cat_raw[[i]] <- sum(mu_cat[[i]]\$False_double_variant_count) / sum(mu_cat[[i]]\$Exposure_m)
    mu_cat_raw[[i]] <- pmax(mu_cat_raw[[i]], 1e-12)
    mu_cat[[i]] <- c(sum(mu_cat[[i]]\$False_double_variant_count) + nu_global * mu_global) / c(sum(mu_cat[[i]]\$Exposure_m) + nu_global)
    mu_cat[[i]] <- pmax(mu_cat[[i]], 1e-12)
  }
  mu_cat_raw <- do.call(c, mu_cat_raw)
  mu_cat <- do.call(c, mu_cat)

  ##### category-specific prior strength nu_cat (ChatGPT function, adapted)
  ## with fallback to global if category is too small
  nu_cat <- rep(NA, 6)
  names(nu_cat) <- names(mu_cat)
  for (i in 1:length(mu_cat)) {
    cat_name <- names(mu_cat)[i]
    idx <- out.summary\$Category == cat_name

    yk <- out.summary\$False_double_variant_count[idx]
    mk <- out.summary\$Exposure_m[idx]
    mu_k <- mu_cat[i]

    opt_k <- try(optimize(f = gp_loglik_nu,
                          interval = c(log(1e-6), log(1e6)),
                          y = yk,
                          m = mk,
                          mu = mu_k),
                 silent = TRUE)

    if (inherits(opt_k, "try-error")) {
      nu_cat[i] <- nu_global
    } else {
      nu_cat[i] <- exp(opt_k\$minimum)
    }
  }

  #### posterior mean for each specific SNV
  out.summary\$mu_category <- rep(NA, nrow(out.summary))
  out.summary\$nu_category <- rep(NA, nrow(out.summary))
  for(i in 1:6){
    out.summary\$mu_category[which(out.summary\$Category == names(mu_cat)[i])] <- mu_cat[i]
    out.summary\$nu_category[which(out.summary\$Category == names(nu_cat)[i])] <- nu_cat[i]
  }

  ##### prior: e_j ~ Gamma(shape = nu_cat * mu_cat, rate = nu_cat)
  out.summary\$prior_shape <- out.summary\$mu_category * out.summary\$nu_category
  out.summary\$prior_rate <- out.summary\$nu_category

  ##### posterior: e_j | data ~ Gamma(shape = prior_shape + F, rate = prior_rate [= nu_cat] + M)
  out.summary\$post_shape <- out.summary\$prior_shape + out.summary\$False_double_variant_count
  out.summary\$post_rate  <- out.summary\$nu_category + out.summary\$Exposure_m

  ##### posterior mean = shrunk EB estimate
  out.summary\$e_EB  <- out.summary\$post_shape / out.summary\$post_rate

  ## 6.) Apply correction factors systematically (use upper 95% percentile to correct, a.k.a. "conservative")
  seq.error.rate[match(out.summary\$SNV, rownames(seq.error.rate)),"1nt false counts/coverage"] <- out.summary\$e_EB
  input.counts.processed[match(out.summary\$SNV, input.counts.processed\$base_mut),"counts_per_cov_corrected"] <- input.counts.processed[match(out.summary\$SNV, input.counts.processed\$base_mut),"counts_per_cov"] - out.summary\$e_EB
  input.counts.processed[match(out.summary\$SNV, input.counts.processed\$base_mut),"counts_corrected"] <- input.counts.processed[match(out.summary\$SNV, input.counts.processed\$base_mut),"counts"] - c(out.summary\$e_EB * input.counts.processed[match(out.summary\$SNV, input.counts.processed\$base_mut),"cov"])

  ## round to nearest integer, do not allow for negative counts
  input.counts.processed[match(out.summary\$SNV, input.counts.processed\$base_mut),"counts_corrected"] <- round(input.counts.processed[match(out.summary\$SNV, input.counts.processed\$base_mut),"counts_corrected"])
  if(any(which(input.counts.processed[,"counts_corrected"] < 0))){
    input.counts.processed[which(input.counts.processed[,"counts_corrected"] < 0),"counts_per_cov_corrected"] <- 0
    input.counts.processed[which(input.counts.processed[,"counts_corrected"] < 0),"counts_corrected"] <- 0
  }

  # Write the processed data
  write.csv(seq.error.rate, file = seq_error_rate_path, row.names = T)
  write.csv(input.counts.processed, file = output_file_path, row.names = F)

}

## -------------------------------------------------------------------------------------------------
## I/O wrappers preserved from the pipeline module (not part of Maxi's estimator logic)
## -------------------------------------------------------------------------------------------------

# propagate corrected counts onto the completed library table (used by e.g. logdiff)
correct_library_completed <- function(completed_path, corrected_filtered_path, output_file_path){
  completed <- read.csv(completed_path)
  corrected <- read.csv(corrected_filtered_path)
  idx <- match(completed\$codon_mut, corrected\$codon_mut)
  has <- !is.na(idx)
  completed\$counts[has] <- corrected\$counts[idx[has]]
  completed\$counts_per_cov[has] <- corrected\$counts_per_cov[idx[has]]
  write.csv(completed, file = output_file_path, row.names = FALSE)
}

# re-make the AA-level heatmap table from corrected counts (canonical column names)
seq_error_correct_counts_for_heatmaps <- function(input_csv_path, aa_seq_file_path, output_csv_path, threshold = 3) {

  raw_counts <- read.table(input_csv_path, sep = ",", header = TRUE, stringsAsFactors = FALSE)

  wt_seq <- readLines(aa_seq_file_path)
  wt_seq <- unlist(strsplit(wt_seq, ""))
  wt_seq[wt_seq == "X"] <- "*"

  # summarise corrected counts for each unique pos_mut
  aggregated_counts_per_cov <- aggregate(
    counts_per_cov_corrected ~ pos_mut,
    data = raw_counts,
    FUN = function(x) sum(x, na.rm = TRUE)
  )
  aggregated_counts <- aggregate(
    counts_corrected ~ pos_mut,
    data = raw_counts,
    FUN = function(x) sum(x, na.rm = TRUE)
  )
  aggregated_data <- merge(aggregated_counts_per_cov, aggregated_counts, by = "pos_mut", all = TRUE)

  # canonical output names so the count heatmaps consume this drop-in
  names(aggregated_data)[names(aggregated_data) == "counts_per_cov_corrected"] <- "total_counts_per_cov"
  names(aggregated_data)[names(aggregated_data) == "counts_corrected"] <- "total_counts"

  aggregated_data\$wt_aa <- sub("(\\\\D)(\\\\d+)(\\\\D)", "\\\\1", aggregated_data\$pos_mut)
  aggregated_data\$position <- as.numeric(sub("(\\\\D)(\\\\d+)(\\\\D)", "\\\\2", aggregated_data\$pos_mut))
  aggregated_data\$mut_aa <- sub("(\\\\D)(\\\\d+)(\\\\D)", "\\\\3", aggregated_data\$pos_mut)
  aggregated_data\$mut_aa[aggregated_data\$mut_aa == "X"] <- "*"

  all_amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                       "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "*")
  all_positions <- seq_along(wt_seq)
  complete_data <- expand.grid(mut_aa = all_amino_acids, position = all_positions, stringsAsFactors = FALSE)

  heatmap_data <- merge(complete_data, aggregated_data, by = c("mut_aa", "position"), all.x = TRUE, sort = FALSE)

  heatmap_data\$total_counts_per_cov[is.na(heatmap_data\$total_counts_per_cov)] <- 0
  heatmap_data\$wt_aa <- wt_seq[heatmap_data\$position]

  low_count <- is.na(heatmap_data\$total_counts) | heatmap_data\$total_counts < threshold
  heatmap_data\$total_counts_per_cov[low_count] <- NA
  heatmap_data\$total_counts[low_count] <- NA

  missing_pos_mut <- is.na(heatmap_data\$pos_mut)
  heatmap_data\$pos_mut[missing_pos_mut] <- paste0(
    heatmap_data\$wt_aa[missing_pos_mut],
    heatmap_data\$position[missing_pos_mut],
    heatmap_data\$mut_aa[missing_pos_mut]
  )

  out.order <- paste0(rep(1:max(heatmap_data\$position), each = 21),
                      rep(c("A", "C", "D", "E", "F", "G", "H",
                            "I", "K", "L","M", "N", "P", "Q",
                            "R", "S", "T", "V", "W", "Y", "*"),
                          max(heatmap_data\$position)))
  heatmap_data <- heatmap_data[match(out.order, paste0(heatmap_data\$position, heatmap_data\$mut_aa)),]
  rownames(heatmap_data) <- 1:nrow(heatmap_data)

  write.csv(heatmap_data, file = output_csv_path, row.names = FALSE)

}

# Defensive coordinate-contract check (general, not GID1A-specific): base_mut positions must be
# 1-based within the supplied reference, and the stated WT base must match the reference base there.
# A mismatch means the reference and the counts are in different frames -> stop rather than
# silently mis-index seq.error.rate.
check_reference_frame <- function(wt_path, processed_path){
  ref <- strsplit(as.character(readDNAStringSet(wt_path)[[1]]), "")[[1]]
  proc <- read.csv(processed_path)
  snv <- proc[grep("[,]", proc[,"base_mut"], invert = TRUE), "base_mut"]
  if(length(snv) == 0) return(invisible(NULL))
  pos <- suppressWarnings(as.integer(str_split_fixed(snv, ":", 2)[,1]))
  wtbase <- substr(str_split_fixed(snv, ":", 2)[,2], 1, 1)
  ok <- !is.na(pos) & pos >= 1 & pos <= length(ref)
  if(any(!ok)){
    stop("Reference/counts coordinate mismatch: ", sum(!ok), " single-nucleotide base_mut positions ",
         "fall outside the reference (length ", length(ref), "). The reference fasta must be the same ",
         "one the reads were aligned to, in the same coordinate frame as base_mut.")
  }
  mism <- ref[pos] != wtbase
  if(any(mism, na.rm = TRUE)){
    stop("Reference/counts coordinate mismatch: ", sum(mism, na.rm = TRUE), " single-nucleotide ",
         "base_mut entries state a wild-type base that disagrees with the reference at that position ",
         "(e.g. ", snv[which(mism)[1]], "). Check that --fasta matches the counts.")
  }
  invisible(NULL)
}

## -------------------------------------------------------------------------------------------------
## Run
## -------------------------------------------------------------------------------------------------
method <- "$method"
codon_window <- $codon_window

# fail loudly if the supplied reference and the counts are in different coordinate frames
check_reference_frame("$wt_fasta", "$filtered")

intermediate <- "variantCounts_intermediate_corrected.csv"
if (method == "eb") {
  seq_error_correct_by_false_doubles_EB(
    wt_path = "$wt_fasta",
    input_count_path_raw = "$raw_counts",
    input_count_path_processed = "$filtered",
    output_file_path = intermediate,
    seq_error_rate_path = "seq_error_rate.csv",
    codon_window = codon_window
  )
} else {
  seq_error_correct_by_false_doubles_MLE(
    wt_path = "$wt_fasta",
    input_count_path_raw = "$raw_counts",
    input_count_path_processed = "$filtered",
    output_file_path = intermediate,
    seq_error_rate_path = "seq_error_rate.csv",
    codon_window = codon_window
  )
}

## corrected counts become the canonical columns so every count-dependent step uses them;
## the originals are preserved as *_raw for the error-correction report
ic <- read.csv(intermediate)
ic\$counts_raw <- ic\$counts
ic\$counts_per_cov_raw <- ic\$counts_per_cov
ic\$counts <- ic\$counts_corrected
ic\$counts_per_cov <- ic\$counts_per_cov_corrected
write.csv(ic, file = "variantCounts_filtered_by_library_error_corrected.csv", row.names = FALSE)

seq_error_correct_counts_for_heatmaps(
  input_csv_path = "variantCounts_filtered_by_library_error_corrected.csv",
  aa_seq_file_path = "$aa_seq",
  output_csv_path = "variantCounts_for_heatmaps_error_corrected.csv",
  threshold = $min_counts
)

correct_library_completed(
  completed_path = "$completed",
  corrected_filtered_path = "variantCounts_filtered_by_library_error_corrected.csv",
  output_file_path = "library_completed_variantCounts_error_corrected.csv"
)

## --- versions ---
r_version <- strsplit(version[['version.string']], ' ')[[1]][3]
if (is.null(r_version)) r_version <- "unknown"
stringr_version <- as.character(packageVersion("stringr"))
biostrings_version <- as.character(packageVersion("Biostrings"))
f <- file("versions.yml", "w")
writeLines(
  c(
    '"${task.process}":',
    paste('    r-base:', r_version),
    paste('    r-stringr:', stringr_version),
    paste('    r-Biostrings:', biostrings_version)
  ),
  f
)
close(f)
