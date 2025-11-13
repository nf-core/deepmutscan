# input: wildtype nucleotide sequence as string or .fa
# output: first nucleotide position of start codon and last nucleotide position of stopp codon (e.g. "7-2145")

# Load necessary package
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}
suppressMessages(library(Biostrings))

# Function to find the nearest stop codon in the reading frame of the first start codon
# Accepts either a sequence string or a path to a reference FASTA file
find_start_stop <- function(input_string, start_codon = "ATG", stop_codons = c("TAA", "TAG", "TGA")) {
  # Determine if the input is a file path or a raw sequence
  if (file.exists(input_string)) {
    # If the input is a valid file path, read the DNA sequence from the FASTA file
    dna <- readDNAStringSet(input_string)
    sequence <- as.character(dna[[1]]) # Extract the sequence as a character string
  } else {
    # If the input is not a file path, treat it as a raw sequence
    sequence <- input_string
  }

  # Check if the sequence is provided or successfully read
  if (is.null(sequence)) {
    stop("No valid sequence or FASTA file provided.")
  }

  # Convert the sequence to uppercase to ensure matching regardless of case
  sequence <- toupper(sequence)

  # Function to find codons in the sequence
  find_codons <- function(sequence, codon) {
    codon_positions <- gregexpr(codon, sequence)[[1]]
    codon_positions[codon_positions > 0] # Return only positive matches
  }

  # Search for the first start codon
  start_positions <- find_codons(sequence, start_codon)

  if (length(start_positions) == 0 || start_positions[1] == -1) {
    return("No start codon found") # Return message if no start codon is found
  }
  first_start <- start_positions[1] # Consider only the first start codon

  # Search for all possible stop codons
  stop_positions_list <- lapply(stop_codons, find_codons, sequence = sequence)
  all_stop_positions <- sort(unlist(stop_positions_list))

  # Filter only stop codons that are in the same reading frame and come after the start codon
  in_frame_stops <- all_stop_positions[(all_stop_positions - first_start) %% 3 == 0 & all_stop_positions > first_start]
  if (length(in_frame_stops) > 0) {
    # Choose the nearest stop codon
    closest_stop <- min(in_frame_stops)
    return(cat(paste(first_start, closest_stop + 2, sep = "-"))) # Output format: start-stop
  } else {
    return("No suitable stop codon found") # If no suitable stop codon is found
  }
}

#find_start_stop("/Users/benjaminwehnert/CRG/DMS_QC/testing_data/MORtn5_reference.fa")


