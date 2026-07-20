#!/usr/bin/env Rscript

# input: wildtype-seq, start&stopp pos.
# output: amino acid sequence within the start-stop frame (.txt)

# Load necessary libraries
suppressMessages(library(Biostrings))

# Define the function
aa_seq <- function(wt_seq_input, pos_range, output_file) {
  # Parse the start and stop positions from the input format "start-stop"
  positions <- unlist(strsplit(pos_range, "-"))
  start_pos <- as.numeric(positions[1])
  stop_pos <- as.numeric(positions[2])

  # Check if the input is a file or a string
  if (file.exists(wt_seq_input)) {
    # If it's a file, read the sequence from the fasta file
    seq_data <- readDNAStringSet(filepath = wt_seq_input)
    wt_seq <- seq_data[[1]]  # Extract the sequence
  } else {
    # Otherwise, treat the input as a sequence string
    wt_seq <- DNAString(wt_seq_input)
  }

  # Extract the sequence between start and stop codons
  coding_seq <- subseq(wt_seq, start = start_pos, end = stop_pos)

  # Translate the coding sequence into an amino acid sequence
  aa_seq <- translate(coding_seq)

  # Write the amino acid sequence to a .txt file
  write(as.character(aa_seq), file = output_file)
}


#####
# run function
#####
aa_seq(
  wt_seq_input = "$wt_seq",
  "$pos_range",
  "aa_seq.txt"
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
    '"${task.process}":',
    paste('    r-base:', r_version),
    paste('    r-Biostrings:', Biostrings_version)
  ),
  f
)
close(f)
