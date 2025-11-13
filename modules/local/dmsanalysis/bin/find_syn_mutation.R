suppressMessages(library(Biostrings))

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
  wt_aa  <- Biostrings::translate(wt_subseq, if.fuzzy.codon = "X")
  wt_chars <- strsplit(wt_seq_chr, "", fixed = TRUE)[[1]]
  
  ## ---- counts ----
  df <- utils::read.delim(counts_merged_tsv, sep = "\t", header = TRUE,
                          stringsAsFactors = FALSE, check.names = FALSE)
  if (!"nt_seq" %in% names(df)) stop("counts_merged_tsv must have a 'nt_seq' column.")
  
  df$nt_seq <- toupper(df$nt_seq)
  keep_len <- nchar(df$nt_seq) == wt_len
  if (!any(keep_len)) stop("No sequences match WT window length (", wt_len, ").")
  if (!all(keep_len)) df <- df[keep_len, , drop = FALSE]
  
  # input columns & mean (works with 1+ replicates)
  input_cols <- grep("^input", names(df), value = TRUE)
  if (length(input_cols) == 0L) stop("No input columns found (expect names starting with 'input').")
  input_mat <- as.data.frame(lapply(df[, input_cols, drop = FALSE], function(x) as.numeric(as.character(x))))
  input_mean <- if (length(input_cols) == 1L) input_mat[[1]] else rowMeans(as.matrix(input_mat), na.rm = TRUE)
  
  ## ---- synonymous filter ----
  var_set <- Biostrings::DNAStringSet(df$nt_seq)
  var_aa  <- Biostrings::translate(var_set, if.fuzzy.codon = "X")
  syn_idx <- which(as.character(var_aa) == as.character(wt_aa))
  if (length(syn_idx) == 0L) stop("No fully-synonymous variants found relative to WT translation.")
  
  # helpers
  mismatch_positions <- function(seq_nt_chars) which(seq_nt_chars != wt_chars) # 1-based positions
  codon_index <- function(pos_vec) floor((pos_vec - 1L) / 3L)                  # 0-based codon bin
  
  # preference 1: exactly 2 mismatches, both within the same codon
  cand_two_one <- Filter(function(i) {
    vchars <- strsplit(df$nt_seq[i], "", fixed = TRUE)[[1]]
    pos <- mismatch_positions(vchars)
    length(pos) == 2L && length(unique(codon_index(pos))) == 1L
  }, syn_idx)
  
  choose_best <- function(idx_vec) idx_vec[ which.max(input_mean[idx_vec]) ]
  
  if (length(cand_two_one) > 0L) {
    best_i <- choose_best(cand_two_one)
    return(as.character(df$nt_seq[best_i]))
  }
  
  # preference 2 (fallback): exactly 1 mismatch (still synonymous)
  cand_one <- Filter(function(i) {
    vchars <- strsplit(df$nt_seq[i], "", fixed = TRUE)[[1]]
    length(mismatch_positions(vchars)) == 1L
  }, syn_idx)
  
  if (length(cand_one) > 0L) {
    best_i <- choose_best(cand_one)
    return(as.character(df$nt_seq[best_i]))
  }
  
  stop("No suitable synonymous variant found: neither 2-in-1-codon nor 1-nt synonymous candidates present.")
}
