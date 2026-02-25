#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Script: Generate Programmed Codon Variants
# Description: Generates all possible programmed codon mutations for a given 
#              wild-type sequence based on a specified mutagenesis strategy.
# Input: 
#   - wt_seq_input: Wild-type sequence (string or path to FASTA file).
#   - start_stop_pos: Target sequence range format "start-stop".
#   - mutagenesis_type: Strategy ('nnk', 'nns', 'max_diff_to_wt', 'custom').
#   - custom_codon_library_path: Path to custom library. Automatically detects 
#     if the file is a global list ("AAA, AAC...") or a position-wise CSV 
#     (requires a "Position" header).
#   - output_file: Desired name/path for the output CSV.
# Output: A CSV file containing all possible mutated codons per position.
# ------------------------------------------------------------------------------

suppressMessages(library(Biostrings))
suppressMessages(library(methods))

generate_possible_variants <- function(wt_seq_input, start_stop_pos, mutagenesis_type,
                                       custom_codon_library_path, output_file) {
  
  # Parse the start and stop positions from the input format "start-stop"
  positions <- unlist(strsplit(start_stop_pos, "-"))
  start_pos <- as.numeric(positions[1])
  stop_pos <- as.numeric(positions[2])
  
  # Load sequence from file or process as a direct string
  if (file.exists(wt_seq_input)) {
    seq_data <- Biostrings::readDNAStringSet(filepath = wt_seq_input)
    wt_seq <- seq_data[[1]]
  } else {
    wt_seq <- Biostrings::DNAString(wt_seq_input)
  }
  
  # Extract the target coding sequence
  coding_seq <- Biostrings::subseq(wt_seq, start = start_pos, end = stop_pos)
  coding_seq <- as.character(coding_seq)
  
  # Predefined codon dictionaries
  nnk_codons <- c('AAG', 'AAT', 'ATG', 'ATT', 'AGG', 'AGT', 'ACG', 'ACT', 'TAG', 'TAT', 'TTG', 'TTT', 'TGG', 'TGT', 'TCG', 'TCT', 'GAG', 'GAT', 'GTG', 'GTT', 'GGG', 'GGT', 'GCG', 'GCT', 'CAG', 'CAT', 'CTG', 'CTT', 'CGG', 'CGT', 'CCG', 'CCT')
  nns_codons <- c('AAG', 'AAC', 'ATG', 'ATC', 'AGG', 'AGC', 'ACG', 'ACC', 'TAG', 'TAC', 'TTG', 'TTC', 'TGG', 'TGC', 'TCG', 'TCC', 'GAG', 'GAC', 'GTG', 'GTC', 'GGG', 'GGC', 'GCG', 'GCC', 'CAG', 'CAC', 'CTG', 'CTC', 'CGG', 'CGC', 'CCG', 'CCC')
  
  # --------------------------------------------------------------------------
  # Custom Library Parsing with Auto-Detection
  # --------------------------------------------------------------------------
  is_position_wise <- FALSE
  position_lookup <- list()
  custom_codons <- NULL
  
  if (mutagenesis_type == "custom") {
    if (!file.exists(custom_codon_library_path) || is.null(custom_codon_library_path)) {
      stop("Custom codons file must be provided and valid when using 'custom' mutagenesis_type.")
    }
    
    # Auto-detect format by inspecting the first line
    first_line <- readLines(custom_codon_library_path, n = 1)
    
    if (grepl("Position", first_line, ignore.case = TRUE)) {
      # Format 1: Position-wise CSV file
      is_position_wise <- TRUE
      
      # Read line-by-line instead of read.csv to avoid strict column matching errors
      lines <- readLines(custom_codon_library_path)
      
      # Loop through lines, skipping the header (index 1)
      for (line in lines[-1]) {
        # Skip empty lines
        if (trimws(line) == "") next 
        
        # Split the line by commas
        parts <- trimws(unlist(strsplit(line, ",")))
        
        # The first part is the position, everything else are the codons
        pos_idx <- parts[1]
        codon_vec <- parts[-1]
        
        # Remove any accidental empty strings (e.g., from trailing commas)
        codon_vec <- codon_vec[codon_vec != ""]
        
        position_lookup[[pos_idx]] <- codon_vec
      }
    } else {
      # Format 2: Global comma-separated list (Legacy compatibility)
      custom_codons <- unlist(strsplit(readLines(custom_codon_library_path), ","))
      custom_codons <- trimws(custom_codons)
    }
  }
  
  # Helper function to split a DNA sequence into nucleotide triplets
  split_into_codons <- function(seq) {
    # Note: Double escaping is required for Perl regular expressions here
    return(strsplit(seq, "(?<=.{3})", perl = TRUE)[[1]])
  }
  
  wt_codons <- split_into_codons(coding_seq)
  
  # Initialize dataframe to store final variant results
  # Note: \$ escaping is maintained for Nextflow compatibility
  result <- data.frame(Codon_Number = integer(), wt_codon = character(), Variant = character(), stringsAsFactors = FALSE)
  
  # Helper function to determine the target codon list per position
  get_codon_list <- function(wt_codon, codon_index) {
    if (mutagenesis_type == "nnk") {
      return(nnk_codons)
    } else if (mutagenesis_type == "nns") {
      return(nns_codons)
    } else if (mutagenesis_type == "max_diff_to_wt") {
      if (substr(wt_codon, 3, 3) == "T") return(nns_codons) else return(nnk_codons)
    } else if (mutagenesis_type == "custom") {
      if (is_position_wise) {
        idx_str <- as.character(codon_index)
        if (!is.null(position_lookup[[idx_str]])) {
          return(position_lookup[[idx_str]])
        } else {
          return(NULL) # Skip positions not explicitly defined in the CSV
        }
      } else {
        return(custom_codons)
      }
    } else {
      stop("Invalid mutagenesis_type. Choose from 'nnk', 'nns', 'max_diff_to_wt', or 'custom'.")
    }
  }
  
  # Iterate over each wild-type codon to assign programmed variants
  for (i in seq_along(wt_codons)) {
    wt_codon <- wt_codons[i]
    codon_list <- get_codon_list(wt_codon, i)
    
    # Skip iteration if no custom codons were assigned to this specific position
    if (is.null(codon_list)) next 
    
    # Filter out the wild-type codon from the mutation list
    possible_variants <- codon_list[codon_list != wt_codon]
    
    for (variant in possible_variants) {
      # Note: \$ escaping is maintained for Nextflow compatibility
      result <- rbind(result, data.frame(Codon_Number = i, wt_codon = wt_codon, Variant = variant, stringsAsFactors = FALSE))
    }
  }
  
  write.csv(result, output_file, row.names = FALSE)
}

# ------------------------------------------------------------------------------
# Main Execution Block (Nextflow variable substitution)
# ------------------------------------------------------------------------------

# Replaces bash if/else logic. If Nextflow omits the optional file, 
# it passes "/NULL", which R translates to an actual NULL object.
custom_lib_arg <- if ("$custom_codon_library" == "/NULL") {
  NULL
} else {
  "$custom_codon_library"
}

generate_possible_variants(
  wt_seq_input              = "$wt_seq",
  start_stop_pos            = "$pos_range",
  mutagenesis_type          = "$mutagenesis_type",
  custom_codon_library_path = "$custom_codon_library",
  output_file               = "possible_mutations.csv"
)

# ------------------------------------------------------------------------------
# Versioning Generation
# ------------------------------------------------------------------------------

r_version <- strsplit(version[['version.string']], ' ')[[1]][3]
biostrings_version <- as.character(packageVersion("Biostrings"))

f <- file("versions.yml", "w")
writeLines(
  c(
    '"${task.process}":',
    paste('    r-base:', r_version),
    paste('    biostrings:', biostrings_version)
  ),
  f
)
close(f)
