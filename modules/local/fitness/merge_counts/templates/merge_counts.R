#!/usr/bin/env Rscript

# --------------------------------------------------------------
# 1. SETUP & INPUTS
# --------------------------------------------------------------

raw_inputs_str  <- "$input_counts"
raw_outputs_str <- "$output_counts"

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
  file = "counts_merged.tsv",
  sep = "\\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

# --------------------------------------------------------------
# 3. VERSIONING
# --------------------------------------------------------------

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
