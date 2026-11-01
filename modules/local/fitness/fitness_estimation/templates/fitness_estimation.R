#!/usr/bin/env Rscript

# Default fitness estimation of one biological sample, in one task: the six steps that used to be six
# processes (MERGE_COUNTS, EXPDESIGN_FITNESS, FIND_SYNONYMOUS_MUTATION, FITNESS_CALCULATION, FITNESS_QC,
# FITNESS_HEATMAP), run in sequence in one R session on the same image. The functions below are those
# scripts' code; the count merging, which was top-level script code, is wrapped in merge_counts().
# Every intermediate file is still written (and published) as before, and the next step reads it back
# exactly as the separate process did.
#
# Inputs:  input_counts / output_counts (per-library DiMSum-format count tables of the sample, sorted
#          by replicate), samplesheet_csv, wt_fasta, pos_range ("start-end", 1-based inclusive),
#          wt_seq (wildtype amino-acid sequence, one line), min_counts, output_pseudocount
# Outputs: counts_merged.tsv, experimentalDesign.tsv, synonymous_wt.txt, fitness_estimation.tsv,
#          fitness_estimation_count_correlation.pdf, fitness_estimation_fitness_correlation.pdf,
#          fitness_heatmap.pdf

suppressPackageStartupMessages({
  library(Biostrings)
  library(methods)
  library(dplyr)
  library(ggplot2)
  library(grid)        # for unit()
})

####################################################################################################
# 1. Merge the per-library counts of the sample (counts_merged.tsv)
####################################################################################################

merge_counts <- function(raw_inputs_str, raw_outputs_str, out_path = "counts_merged.tsv") {

  parse_paths <- function(x) {
    if (x == "" || x == "[]") return(character(0))
    # Split at spaces
    parts <- strsplit(x, "\\\\s+")[[1]]
    # Remove empty strings
    parts[parts != ""]
  }

  input_paths  <- parse_paths(raw_inputs_str)
  output_paths <- parse_paths(raw_outputs_str)

  message("Found ", length(input_paths),  " input files.")
  message("Found ", length(output_paths), " output files.")

  # --------------------------------------------------------------
  # 2. MERGE LOGIC
  # --------------------------------------------------------------

  # Helper to read a 2-col TSV without header
  read_counts <- function(fp) {
    utils::read.table(
      fp, header = FALSE, sep = "\\t", quote = "",
      col.names = c("nt_seq", "count"),
      colClasses = c("character", "numeric"),
      comment.char = "", check.names = FALSE
    )
  }

  # Read all inputs / outputs
  input_list  <- lapply(input_paths,  read_counts)
  output_list <- lapply(output_paths, read_counts)

  # Collect universe of sequences
  all_seqs <- unique(c(
    unlist(lapply(input_list,  function(x) x\$nt_seq)),
    unlist(lapply(output_list, function(x) x\$nt_seq))
  ))

  # Pre-allocate output frame
  n_in  <- length(input_list)
  n_out <- length(output_list)

  col_names <- c(
    "nt_seq",
    if (n_in  > 0) paste0("input",  seq_len(n_in))  else character(0),
    if (n_out > 0) paste0("output", seq_len(n_out)) else character(0)
  )

  # Initialize dataframe with 0 counts
  out <- data.frame(
    nt_seq = all_seqs,
    matrix(0, nrow = length(all_seqs), ncol = n_in + n_out),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  names(out) <- col_names

  # Fill inputs
  if (n_in > 0) {
    for (i in seq_len(n_in)) {
      df <- input_list[[i]]
      # Match sequences: WICHTIG \$ escaping
      idx <- match(df\$nt_seq, out\$nt_seq)
      # Assign counts: WICHTIG \$ escaping
      out[idx, paste0("input", i)] <- df\$count
    }
  }

  # Fill outputs
  if (n_out > 0) {
    for (j in seq_len(n_out)) {
      df <- output_list[[j]]
      # WICHTIG \$ escaping
      idx <- match(df\$nt_seq, out\$nt_seq)
      out[idx, paste0("output", j)] <- df\$count
    }
  }

  # Write Output
  utils::write.table(
    out,
    file = out_path,
    sep = "\\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
}

####################################################################################################
# 2. DiMSum experimental design from the samplesheet (experimentalDesign.tsv)
####################################################################################################

# Make a DiMSum experimental design from a deepmutscan samplesheet.
# - samplesheet_csv: path to CSV with columns sample,type,replicate,file1,file2
# - out_path: where to write the TSV (default "experimentalDesign.tsv")
# Returns: the experimental design as a data.frame
make_dimsum_experimental_design <- function(samplesheet_csv, out_path = "experimentalDesign.tsv") {
  # ---- read & normalize ----
  ss <- read.csv(samplesheet_csv, stringsAsFactors = FALSE, check.names = FALSE)
  names(ss) <- tolower(names(ss))

  # tolerate missing file2 column (single-end)
  if (!"file2" %in% names(ss)) ss\$file2 <- ""

  required <- c("sample", "type", "replicate", "file1", "file2")
  missing  <- setdiff(required, names(ss))
  if (length(missing) > 0) stop("Samplesheet missing columns: ", paste(missing, collapse = ", "))

  # coerce types
  ss\$replicate <- as.integer(ss\$replicate)

  # ---- keep only the samples DiMSum scores ----
  # A samplesheet legitimately carries other row types: `wildtype` (required by
  # --error_correction wildtype) and `quality`. They are real samples, but they are not part of the
  # selection design - they feed error correction and QC - and MERGE_COUNTS drops them too. Leaving
  # them in gave them an empty `selection_id`, which DiMSum reads back as NA and then dies in
  # dimsum__check_experiment_design on `sum(exp_design[,"selection_id"] < 0)` with the opaque
  # "missing value where TRUE/FALSE needed".
  ss <- ss[ss\$type %in% c("input", "output"), , drop = FALSE]
  if (nrow(ss) == 0) {
    stop("Samplesheet has no 'input'/'output' rows - DiMSum needs both to compute fitness.")
  }

  # ---- order rows: by sample (if multiple), input before output, then replicate ----
  # Sorting happens before naming, because the names are derived from the row order below.
  multi_base <- length(unique(ss\$sample)) > 1
  type_rank  <- match(ss\$type, c("input", "output"))
  ord <- if (multi_base) {
    order(ss\$sample, type_rank, ss\$replicate)
  } else {
    order(type_rank, ss\$replicate)
  }
  ss <- ss[ord, , drop = FALSE]

  # ---- sample_name: must equal the count column MERGE_COUNTS writes for this row ----
  # merge_counts.R names its columns positionally - paste0("input", seq_len(n_in)) over the files
  # sorted by replicate - so the name has to come from each row's POSITION within its type, not from
  # the samplesheet's `replicate` value. The two only coincide when the replicates happen to be a
  # contiguous 1..N; replicates 1,2,4 would otherwise produce "input4" for a column named "input3"
  # and DiMSum could not match them. DiMSum also requires sample_name to be alphanumeric and to
  # start with a letter, so the parts are joined without a separator.
  idx <- stats::ave(seq_len(nrow(ss)), ss\$sample, ss\$type, FUN = seq_along)
  # NOTE: for multiple biological samples these names still will not match the counts file, because
  # one design is built from the whole samplesheet and broadcast to every per-sample DiMSum run
  # (see subworkflows/local/calculate_fitness/main.nf). Multi-sample + DiMSum is a separate,
  # known-broken case; single-sample runs are correct.
  sample_name <- if (multi_base) paste0(ss\$sample, ss\$type, idx) else paste0(ss\$type, idx)

  # ---- build DiMSum columns ----
  experiment_replicate <- ss\$replicate           # pairs an input with an output; need not be 1..N
  selection_id <- ifelse(ss\$type == "input", 0L, 1L)
  # assume one selection batch; DiMSum requires this to be blank for input samples
  selection_replicate <- ifelse(ss\$type == "output", 1L, NA_integer_)
  # assume one technical batch
  technical_replicate <- rep(1L, nrow(ss))

  pair1 <- basename(ss\$file1)
  # keep empty string for single-end / missing file2
  pair2 <- ifelse(is.na(ss\$file2) | ss\$file2 == "", "", basename(ss\$file2))

  ed <- data.frame(
    sample_name          = sample_name,
    experiment_replicate = experiment_replicate,
    selection_id         = selection_id,
    selection_replicate  = selection_replicate,
    technical_replicate  = technical_replicate,
    pair1                = pair1,
    pair2                = pair2,
    stringsAsFactors     = FALSE
  )
  rownames(ed) <- NULL

  # ---- write & return ----
  write.table(ed, file = out_path, sep = "\\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")
  return(ed)
}

####################################################################################################
# 3. Synonymous wildtype proxy for the normalisation (synonymous_wt.txt)
####################################################################################################

# Pick a synonymous "WT substitute" for DiMSum normalization using a fixed coding window.
# Inputs:
#   wt_fasta          : path to FASTA (single WT sequence)
#   counts_merged_tsv : path to merged counts (columns: nt_seq, input1..N, output1..M)
#   pos_range         : "start-end" (1-based, inclusive), e.g. "352-1383"
# Returns:
#   character scalar: chosen nt sequence
#
# Preference:
#   1) AA-identical (fully synonymous) AND exactly 2 nt mismatches vs WT, both within ONE codon.
#   2) If none, AA-identical AND exactly 1 nt mismatch vs WT.
#   3) If more than one, pick highest mean of input counts.
#   4) If still none, stop with an error.

pick_synonymous_wt_from_range <- function(wt_fasta, counts_merged_tsv, pos_range) {
  ## ---- parse range ----
  pr <- strsplit(as.character(pos_range), "-", fixed = TRUE)[[1]]
  if (length(pr) != 2L) stop("pos_range must be 'start-end', got: ", pos_range)
  start_pos <- as.integer(pr[1]); end_pos <- as.integer(pr[2])
  if (is.na(start_pos) || is.na(end_pos) || start_pos < 1L || end_pos < start_pos)
    stop("Invalid pos_range: ", pos_range)

  ## ---- WT window ----
  wt_set <- Biostrings::readDNAStringSet(wt_fasta)
  if (length(wt_set) != 1L) stop("WT FASTA must contain exactly one sequence.")
  wt_subseq <- Biostrings::subseq(wt_set[[1]], start = start_pos, end = end_pos)
  wt_seq_chr <- as.character(wt_subseq)
  wt_len <- nchar(wt_seq_chr)
  if ((wt_len %% 3) != 0) stop("Provided window length is not divisible by 3: ", wt_len)
  wt_aa  <- Biostrings::translate(wt_subseq, if.fuzzy.codon = "X", no.init.codon = TRUE)
  wt_chars <- strsplit(wt_seq_chr, "", fixed = TRUE)[[1]]

  ## ---- counts ----
  df <- utils::read.delim(counts_merged_tsv, sep = "\\t", header = TRUE,
                          stringsAsFactors = FALSE, check.names = FALSE)
  if (!"nt_seq" %in% names(df)) stop("counts_merged_tsv must have a 'nt_seq' column.")

  df\$nt_seq <- toupper(df\$nt_seq)
  keep_len <- nchar(df\$nt_seq) == wt_len
  if (!any(keep_len)) stop("No sequences match WT window length (", wt_len, ").")
  if (!all(keep_len)) df <- df[keep_len, , drop = FALSE]

  # input columns & mean (works with 1+ replicates)
  input_cols <- grep("^input", names(df), value = TRUE)
  if (length(input_cols) == 0L) stop("No input columns found (expect names starting with 'input').")
  input_mat <- as.data.frame(lapply(df[, input_cols, drop = FALSE], function(x) as.numeric(as.character(x))))
  input_mean <- if (length(input_cols) == 1L) input_mat[[1]] else rowMeans(as.matrix(input_mat), na.rm = TRUE)

  ## ---- synonymous filter ----
  var_set <- Biostrings::DNAStringSet(df\$nt_seq)
  var_aa  <- Biostrings::translate(var_set, if.fuzzy.codon = "X", no.init.codon = TRUE)
  syn_idx <- which(as.character(var_aa) == as.character(wt_aa))
  if (length(syn_idx) == 0L) stop("No fully-synonymous variants found relative to WT translation.")

  # helpers
  mismatch_positions <- function(seq_nt_chars) which(seq_nt_chars != wt_chars) # 1-based positions
  codon_index <- function(pos_vec) floor((pos_vec - 1L) / 3L)                  # 0-based codon bin

  # preference 1: exactly 2 mismatches, both within the same codon
  cand_two_one <- Filter(function(i) {
    vchars <- strsplit(df\$nt_seq[i], "", fixed = TRUE)[[1]]
    pos <- mismatch_positions(vchars)
    length(pos) == 2L && length(unique(codon_index(pos))) == 1L
  }, syn_idx)

  choose_best <- function(idx_vec) idx_vec[ which.max(input_mean[idx_vec]) ]

  if (length(cand_two_one) > 0L) {
    best_i <- choose_best(cand_two_one)
    return(as.character(df\$nt_seq[best_i]))
  }

  # preference 2 (fallback): exactly 1 mismatch (still synonymous)
  cand_one <- Filter(function(i) {
    vchars <- strsplit(df\$nt_seq[i], "", fixed = TRUE)[[1]]
    length(mismatch_positions(vchars)) == 1L
  }, syn_idx)

  if (length(cand_one) > 0L) {
    best_i <- choose_best(cand_one)
    return(as.character(df\$nt_seq[best_i]))
  }

  stop("No suitable synonymous variant found: neither 2-in-1-codon nor 1-nt synonymous candidates present.")
}

####################################################################################################
# 4. Default fitness estimation (fitness_estimation.tsv)
####################################################################################################

## default fitness estimation for nf-core/deepmutscan
## 18.02.2026
## maximilian.stammnitz@crg.eu

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
  merged.counts <- cbind("aa_seq" = as.character(translate(DNAStringSet(merged.counts\$nt_seq), no.init.codon = TRUE)), merged.counts)
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

# keep nucleotide variants with at least `min_counts` reads in every input replicate, before
# synonymous codons are merged - the same filter as DiMSum's '--fitnessMinInputCountAll'
# (the wildtype is always kept, as it anchors every fitness estimate)
filter_low_input_counts <- function(merged.counts, wt.seq, min_counts) {
  input.cols <- grep("^input[0-9]+\$", colnames(merged.counts), value = TRUE)
  keep <- apply(merged.counts[, input.cols, drop = FALSE] >= min_counts, 1, all)
  keep[is.na(keep)] <- FALSE
  keep[merged.counts\$nt_seq == as.character(wt.seq)] <- TRUE
  merged.counts[keep, , drop = FALSE]
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
calc_raw_fitness <- function(merged.counts, exp.design, output_pseudocount) {
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

    ### dropouts (zero output count, non-zero input count): add the output pseudocount,
    ### or leave them without a fitness estimate in this replicate if it is 0
    tmp.dropout <- which(tmp.output.counts == 0 & tmp.input.counts != 0)
    if (output_pseudocount > 0) {
      tmp.output.counts[tmp.dropout] <- output_pseudocount
    } else {
      tmp.output.counts[tmp.dropout] <- NA
    }

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
    rm(tmp.fitness, tmp.wt.log.ratio, tmp.output.counts, tmp.input.counts, tmp.dropout)
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
#' @param min_counts Minimum number of reads in every input replicate for a fitness estimate
#' @param output_pseudocount Pseudocount for output dropouts (0 = no estimate for dropouts)
#'
#' @return Invisibly returns the final data.frame; writes the output to output_path.
run_fitness_estimation <- function(counts_path,
                                   design_path,
                                   wt_seq_path,
                                   output_path,
                                   min_counts = 10,
                                   output_pseudocount = 1) {
  ## 1. Import key files ##
  #########################

  merged.counts <- read.table(counts_path, sep = "\t", header = TRUE, check.names = FALSE)
  exp.design   <- read.table(design_path, sep = "\t", header = TRUE, check.names = FALSE)
  wt.seq       <- DNAString(as.character(read.table(wt_seq_path)))
  wt.seq.aa    <- translate(wt.seq, no.init.codon = TRUE)

  ## 2. Pre-processing the count table ##
  #######################################

  ## keep variants with at least `min_counts` reads in every input replicate
  merged.counts <- filter_low_input_counts(merged.counts, wt.seq, min_counts)

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
  fitness_res <- calc_raw_fitness(merged.counts, exp.design, output_pseudocount)
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


####################################################################################################
# 5. Fitness QC plots (fitness_estimation_count_correlation.pdf, fitness_estimation_fitness_correlation.pdf)
####################################################################################################

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
    ## Axis range from the data instead of a hard-coded c(-3, 1), which silently cropped
    ## strongly deleterious variants out of the replicate-correlation panels. Fitness is
    ## rescaled so synonymous = 0 and stop = -1, so the usual range is roughly [-1.5, 0.5],
    ## but nothing bounds it: padding the observed range keeps every point visible.
    fitness.range <- range(unlist(merged.counts.fitness[, fitness.repl]), finite = TRUE)
    if (!all(is.finite(fitness.range))) {
      fitness.range <- c(-3, 1)
    }
    fitness.pad <- max(0.1, diff(fitness.range) * 0.05)
    fitness.lim <- c(fitness.range[1] - fitness.pad, fitness.range[2] + fitness.pad)

    pairs(merged.counts.fitness[, fitness.repl],
          lower.panel = panel_xy_abline_fitness,
          upper.panel = panel_cor_fitness,
          cex.text = 2,
          xlim = fitness.lim,
          ylim = fitness.lim)
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

####################################################################################################
# 6. Fitness heatmaps (fitness_heatmap.pdf)
####################################################################################################

# ---------- helper functions ----------
find_col <- function(df, candidates) {
  norm <- function(x) gsub("[^a-z0-9]+", "_", tolower(x))
  nms <- colnames(df); nn <- norm(nms)
  for (cand in candidates) {
    hit <- which(nn == norm(cand))
    if (length(hit) == 1) return(nms[hit])
  }
  stop(sprintf("Could not find any of columns: %s", paste(candidates, collapse = ", ")))
}

get_rescaled_cols <- function(df) {
  nms  <- colnames(df)
  hits <- grep("^rescaled[_ ]?fitness", nms, ignore.case = TRUE, value = TRUE)
  if (!length(hits)) stop("No 'rescaled_fitness' columns found.")
  idx <- suppressWarnings(as.integer(gsub(".*?([0-9]+)\$", "\\\\1", hits))) 	# additional backslashs to make it groovy readable
  hits[order(is.na(idx), idx, hits)]
}

# Find "mean fitness" column (if present)
find_mean_col <- function(df) {
  nms <- colnames(df)
  key <- tolower(gsub("[^a-z0-9]+", "_", nms))
  hit <- which(key == "mean_fitness")
  if (length(hit) == 1) nms[hit] else NULL
}

# NEW: read WT amino acid sequence from .txt file (single line)
read_wt_seq_aa_txt <- function(path) {
  if (is.null(path)) stop("wt_seq_aa_txt_path must be provided.")
  x <- readLines(path, warn = FALSE)
  x <- x[nzchar(x)]
  if (!length(x)) stop("WT AA TXT is empty.")
  aa <- toupper(gsub("\\\\s+", "", x[which.max(nchar(x))]))
  # Keep only valid amino acid letters (including stop '*')
  aa <- gsub("[^ACDEFGHIKLMNPQRSTVWY*]", "", aa)
  if (!nchar(aa)) stop("WT AA TXT contains no valid AA letters.")
  aa
}

# Build full AA×position grid for positions 1..wt_len (from WT sequence),
# join fitness values; any missing combos stay NA (-> grey).
# Positions > wt_len are considered padded and will be white.
build_heatmap_long <- function(df,
                               wt_aa_col,
                               pos_col,
                               mut_aa_col,
                               fitness_col,
                               positions_per_row = 75,
                               wt_seq_aa,
                               fill_missing_as_zero = FALSE) {

  # authoritative WT length from provided sequence
  letters <- strsplit(wt_seq_aa, "", fixed = TRUE)[[1]]
  wt_len  <- length(letters)

  # normalize data
  df0 <- df %>%
    transmute(
      position = suppressWarnings(as.numeric(.data[[pos_col]])),
      wt_aa_in = .data[[wt_aa_col]],
      mut_aa   = .data[[mut_aa_col]],
      fitness  = suppressWarnings(as.numeric(.data[[fitness_col]]))
    ) %>%
    filter(is.finite(position))

  # drop any rows that claim positions beyond WT length
  if (nrow(df0) && any(df0\$position > wt_len, na.rm = TRUE)) {
    dropped <- sum(df0\$position > wt_len, na.rm = TRUE)
    warning(sprintf("Dropping %d row(s) with position > WT length (%d).", dropped, wt_len))
    df0 <- df0 %>% filter(position <= wt_len)
  }

  # pad to next multiple of 75 (by rows)
  rem       <- wt_len %% positions_per_row
  pad_need  <- if (rem == 0) 0 else positions_per_row - rem
  max_paded <- wt_len + pad_need

  # full grid: positions 1..max_paded (so the tail exists), AA set of 21
  all_positions   <- seq_len(max_paded)
  all_amino_acids <- c("A","C","D","E","F","G","H","I","K","L",
                       "M","N","P","Q","R","S","T","V","W","Y","*")

  grid_df <- expand.grid(position = all_positions,
                         mut_aa   = all_amino_acids,
                         KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) %>%
    mutate(is_padded = position > wt_len)

  # join fitness only for real positions (<= wt_len)
  fit_df <- df0 %>% select(position, mut_aa, fitness)
  d <- grid_df %>%
    left_join(fit_df, by = c("position","mut_aa"))

  if (fill_missing_as_zero) {
    d\$fitness[is.na(d\$fitness) & d\$position <= wt_len] <- 0
  }

  # authoritative WT AA per real position; tail gets placeholder 'Y'
  wt_map <- tibble(position = seq_len(wt_len), wt_aa = letters)
  d <- d %>%
    left_join(wt_map, by = "position") %>%
    mutate(wt_aa = ifelse(is.na(wt_aa) & position > wt_len, "Y", wt_aa))

  # layout fields
  d <- d %>%
    mutate(
      row_group = ((position - 1) %/% positions_per_row) + 1,
      wt_aa_pos = paste0(wt_aa, position),
      wt_aa_pos = factor(wt_aa_pos, levels = unique(wt_aa_pos)),
      synonymous = mut_aa == wt_aa
    )

  # IMPORTANT: use WT length as the true end of the protein
  d\$max_pos <- wt_len
  d
}

syn_segments <- function(d, positions_per_row = 75) {
  amino_order <- rev(c("G", "A", "V", "L", "M", "I", "F",
                         "Y", "W", "K", "R", "H", "D", "E",
                         "S", "T", "C", "N", "Q", "P", "*"))
  d %>%
    mutate(
      mut_aa = factor(mut_aa, levels = amino_order),
      x = as.numeric(factor(wt_aa_pos, levels = levels(wt_aa_pos))) -
          ((row_group - 1) * positions_per_row),
      y = as.numeric(factor(mut_aa, levels = amino_order))
    ) %>%
    filter(synonymous, position <= max_pos)
}

# Draw one solid white rectangle per row group covering the padded tail region
white_tail_rects <- function(d, positions_per_row = 75) {
  wt_len <- unique(d\$max_pos)[1]
  if (!is.finite(wt_len)) return(dplyr::tibble()[0,])

  # if perfectly divisible by 75, there is no tail to cover
  if (wt_len %% positions_per_row == 0) return(dplyr::tibble()[0,])

  # which facet (row group) contains the last real position?
  last_group     <- ((wt_len - 1) %/% positions_per_row) + 1
  last_local_idx <- ((wt_len - 1) %% positions_per_row) + 1  # 1..75 within the facet

  tibble::tibble(
    row_group = last_group,
    xmin = last_local_idx + 0.5 - 0.025,     # tiny epsilon to avoid hairlines
    xmax = positions_per_row + 0.5 + 0.025,
    ymin = 0.5  - 0.025,
    ymax = 21.5 + 0.025
  )
}

plot_heatmap <- function(d, title_text, positions_per_row = 75) {
  amino_order <- rev(c("G", "A", "V", "L", "M", "I", "F",
                         "Y", "W", "K", "R", "H", "D", "E",
                         "S", "T", "C", "N", "Q", "P", "*"))
  d <- d %>% mutate(mut_aa = factor(mut_aa, levels = amino_order))

  min_f <- suppressWarnings(min(d\$fitness, na.rm = TRUE)); if (!is.finite(min_f)) min_f <- 0
  max_f <- suppressWarnings(max(d\$fitness, na.rm = TRUE)); if (!is.finite(max_f)) max_f <- 0
  max_orig_pos <- unique(d\$max_pos)[1]

  syn  <- syn_segments(d, positions_per_row)
  rect <- white_tail_rects(d, positions_per_row)

  ggplot(d, aes(x = wt_aa_pos, y = mut_aa, fill = fitness)) +
    scale_fill_gradientn(
      colours = c("#D73027", "#F0F0F0", "#4575B4"),
      values  = if ((abs(min_f) + max_f) > 0) c(0, abs(min_f)/(abs(min_f)+max_f), 1) else c(0, 0.5, 1),
      na.value = "grey35",
      limits   = c(min_f, max_f)
    ) +
    scale_x_discrete(
      labels = function(x) {
        num <- suppressWarnings(as.numeric(gsub("[^0-9]", "", x)))
        ifelse(num > max_orig_pos, " ", x)
      },
      expand = expansion(mult = c(0, 0))  # no extra margin area
    ) +
    geom_tile() +
    # Solid white block covering the tail (no pattern / no seams)
    { if (nrow(rect)) geom_rect(data = rect, inherit.aes = FALSE,
                                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                                fill = "white", color = NA) } +
    geom_segment(
      data = syn,
      aes(x = x - 0.485, xend = x + 0.485, y = y - 0.485, yend = y + 0.485),
      linewidth = 0.2, inherit.aes = FALSE, color = "grey10"
    ) +
    theme_minimal() +
    labs(title = title_text, x = "Wild-type amino acid", y = "Mutant amino acid", fill = "Fitness") +
    theme(
      plot.title   = element_text(size = 16, face = "bold"),
      axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
      axis.text.y  = element_text(size = 10),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
      panel.spacing = grid::unit(0.2, "lines")
    ) +
    facet_wrap(~ row_group, scales = "free_x", ncol = 1)
}

# ---------- main callable ----------
# fitness_table_path : path to fitness_estimation.tsv
# wt_seq_aa_txt_path : path to TXT file containing WT AA sequence (one line)
# output_pdf_path    : output PDF (default "fitness_heatmap.pdf")
# positions_per_row  : default 75
run_fitness_rescaled_heatmaps <- function(fitness_table_path,
                                          wt_seq_aa_txt_path,
                                          output_pdf_path = "fitness_heatmap.pdf",
                                          positions_per_row = 75) {

  df <- read.table(
    fitness_table_path, sep = "\t", header = TRUE,
    check.names = FALSE, quote = "", comment.char = ""
  )

  wt_aa_col     <- find_col(df, c("wt aa", "wt_aa", "wt"))
  pos_col       <- find_col(df, c("pos", "position"))
  mut_aa_col    <- find_col(df, c("mut aa", "mut_aa", "aa"))
  rescaled_cols <- get_rescaled_cols(df)

  wt_seq_aa <- read_wt_seq_aa_txt(wt_seq_aa_txt_path)

  plots <- list()

  ## 1) Mean first – use existing "mean fitness" column if available
  mean_col <- find_mean_col(df)
  if (is.null(mean_col)) {
    df\$`rescaled_fitness_mean` <- if (length(rescaled_cols) == 1) df[[rescaled_cols[1]]] else rowMeans(df[, rescaled_cols], na.rm = TRUE)
    mean_col <- "rescaled_fitness_mean"
  }
  long_df_mean <- build_heatmap_long(df, wt_aa_col, pos_col, mut_aa_col, mean_col,
                                     positions_per_row, wt_seq_aa)
  plots[[length(plots) + 1]] <- list(
    title = sprintf("Fitness — mean of %d replicate(s)", length(rescaled_cols)),
    data  = long_df_mean
  )

  ## 2) Then individual replicates
  for (i in seq_along(rescaled_cols)) {
    col <- rescaled_cols[i]
    long_df <- build_heatmap_long(df, wt_aa_col, pos_col, mut_aa_col, col,
                                  positions_per_row, wt_seq_aa)
    plots[[length(plots) + 1]] <- list(
      title = sprintf("Fitness — rep%d", i),
      data  = long_df
    )
  }

  # Device height: (#row groups × 4)
  page_heights <- vapply(plots, function(p) max(p\$data\$row_group, na.rm = TRUE), numeric(1))
  device_height <- max(4, as.numeric(page_heights) * 4, na.rm = TRUE)

  grDevices::pdf(output_pdf_path, width = 16, height = device_height)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  for (p in plots) print(plot_heatmap(p\$data, p\$title, positions_per_row))
  invisible(TRUE)
}

#####
# run functions
#####
merge_counts(
  raw_inputs_str  = "$input_counts",
  raw_outputs_str = "$output_counts",
  out_path        = "counts_merged.tsv"
)

make_dimsum_experimental_design(
  samplesheet_csv = "$samplesheet_csv",
  out_path = "experimentalDesign.tsv"
)

seq <- pick_synonymous_wt_from_range(
	wt_fasta = "$wt_fasta",
	counts_merged_tsv = "counts_merged.tsv",
	pos_range = "$pos_range"
	)
write(seq, file='synonymous_wt.txt')

run_fitness_estimation(
  counts_path = "counts_merged.tsv",
  design_path = "experimentalDesign.tsv",
  wt_seq_path = "synonymous_wt.txt",
  output_path = "fitness_estimation.tsv",
  min_counts = as.numeric("$min_counts"),
  output_pseudocount = as.numeric("$output_pseudocount")
)

run_fitness_plots(
  fitness_table_path = "fitness_estimation.tsv",
  out_counts_corr_pdf = "fitness_estimation_count_correlation.pdf",
  out_fitness_corr_pdf = "fitness_estimation_fitness_correlation.pdf"
)

run_fitness_rescaled_heatmaps(
  fitness_table_path = "fitness_estimation.tsv",
  wt_seq_aa_txt_path = "$wt_seq",
  output_pdf_path    = "fitness_heatmap.pdf"
)

####
# create versions.yml
####
r_version <- strsplit(version[['version.string']], ' ')[[1]][3]
Biostrings_version <- as.character(packageVersion("Biostrings"))
dplyr_version <- as.character(packageVersion("dplyr"))
ggplot2_version <- as.character(packageVersion("ggplot2"))
methods_version <- as.character(packageVersion("methods"))
grid_version <- as.character(packageVersion("grid"))

if (is.null(r_version)) r_version <- "unknown"
if (length(Biostrings_version) == 0) Biostrings_version <- "unknown"
if (length(dplyr_version) == 0) dplyr_version <- "unknown"
if (length(ggplot2_version) == 0) ggplot2_version <- "unknown"
if (length(methods_version) == 0) methods_version <- "unknown"
if (length(grid_version) == 0) grid_version <- "unknown"

f <- file("versions.yml", "w")
writeLines(
  c(
    '"${task.process}":',
    paste('    r-base:', r_version),
    paste('    r-Biostrings:', Biostrings_version),
    paste('    r-dplyr:', dplyr_version),
    paste('    r-ggplot2:', ggplot2_version),
    paste('    r-methods:', methods_version),
    paste('    r-grid:', grid_version)
  ),
  f
)
close(f)
