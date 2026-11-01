#!/usr/bin/env Rscript

# The two reference tables of a run, in one task: the wildtype amino-acid sequence of the reading
# frame (aa_seq.txt) and the table of all programmed codon variants (possible_mutations.csv). They
# used to be two processes (DMSANALYSIS_AASEQ and DMSANALYSIS_POSSIBLE_MUTATIONS) on the same FASTA,
# reading frame and R image; the two functions below are those scripts, with the reference read and
# translated in base R (the Biostrings calls they made, with the same results).
#
# Inputs:  wt_seq (reference FASTA), pos_range ("start-stop", 1-based inclusive), mutagenesis_type,
#          custom_codon_library (a file, or empty when --custom_codon_library is not given)
# Outputs: aa_seq.txt, possible_mutations.csv

####################################################################################################
# Base-R replacements for Biostrings (readDNAStringSet, subseq, translate), with Biostrings' results
####################################################################################################

# The standard genetic code (Biostrings::GENETIC_CODE); stop codons are "*"
GENETIC_CODE_STANDARD <- c(
  TTT = "F", TTC = "F", TTA = "L", TTG = "L", CTT = "L", CTC = "L", CTA = "L", CTG = "L",
  ATT = "I", ATC = "I", ATA = "I", ATG = "M", GTT = "V", GTC = "V", GTA = "V", GTG = "V",
  TCT = "S", TCC = "S", TCA = "S", TCG = "S", CCT = "P", CCC = "P", CCA = "P", CCG = "P",
  ACT = "T", ACC = "T", ACA = "T", ACG = "T", GCT = "A", GCC = "A", GCA = "A", GCG = "A",
  TAT = "Y", TAC = "Y", TAA = "*", TAG = "*", CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
  AAT = "N", AAC = "N", AAA = "K", AAG = "K", GAT = "D", GAC = "D", GAA = "E", GAG = "E",
  TGT = "C", TGC = "C", TGA = "*", TGG = "W", CGT = "R", CGC = "R", CGA = "R", CGG = "R",
  AGT = "S", AGC = "S", AGA = "R", AGG = "R", GGT = "G", GGC = "G", GGA = "G", GGG = "G"
)

# Sequences of a FASTA file, as readDNAStringSet() reads them: one string per record, sequence lines
# joined, upper case, letters outside the DNA alphabet (IUPAC codes, "-", "+", ".") dropped
read_fasta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  starts <- grep("^>", lines)
  if (length(starts) == 0) stop("no FASTA record in ", path)
  ends <- c(starts[-1] - 1L, length(lines))
  seqs <- mapply(function(s, e) if (e > s) paste(lines[(s + 1L):e], collapse = "") else "", starts, ends)
  gsub("[^ACGTMRWSYKVHDBN.+-]", "", toupper(seqs))
}

# subseq(): stops on a range outside the sequence, as Biostrings does
subseq_checked <- function(x, start, end) {
  if (is.na(start) || is.na(end) || start < 1 || end > nchar(x) || end < start - 1)
    stop("Invalid sequence coordinates: ", start, "-", end, " in a sequence of length ", nchar(x))
  substr(x, start, end)
}

# translate(x, no.init.codon = TRUE): the codons of each sequence by the standard code, trailing
# bases beyond the last full codon ignored. A codon with a letter outside ACGT stops the run, as
# Biostrings does (fuzzy = "error"), or becomes "X" (fuzzy = "X", Biostrings' if.fuzzy.codon = "X").
translate_dna <- function(x, fuzzy = c("error", "X")) {
  fuzzy <- match.arg(fuzzy)
  x <- toupper(x)
  n <- nchar(x) %/% 3L
  out <- character(length(x))
  for (len in unique(n)) {
    i <- which(n == len)
    if (len == 0L) next
    aa <- lapply(seq_len(len), function(k) GENETIC_CODE_STANDARD[substr(x[i], 3L * k - 2L, 3L * k)])
    bad <- vapply(aa, anyNA, logical(1))
    if (any(bad)) {
      if (fuzzy == "error") stop("not a base in codon ", which(bad)[1], " of a sequence to translate")
      aa <- lapply(aa, function(a) { a[is.na(a)] <- "X"; a })
    }
    out[i] <- do.call(paste0, unname(aa))
  }
  out
}

####################################################################################################
# 1. Wildtype amino-acid sequence (aa_seq.txt)
####################################################################################################
# input: wildtype-seq, start&stopp pos.
# output: amino acid sequence within the start-stop frame (.txt)

# Define the function
aa_seq <- function(wt_seq_input, pos_range, output_file) {
  # Parse the start and stop positions from the input format "start-stop"
  positions <- unlist(strsplit(pos_range, "-"))
  start_pos <- as.numeric(positions[1])
  stop_pos <- as.numeric(positions[2])

  # Check if the input is a file or a string
  if (file.exists(wt_seq_input)) {
    # If it's a file, read the sequence from the fasta file
    wt_seq <- read_fasta(wt_seq_input)[1]  # the first sequence
  } else {
    # Otherwise, treat the input as a sequence string
    wt_seq <- toupper(wt_seq_input)
  }

  # Extract the sequence between start and stop codons
  coding_seq <- subseq_checked(wt_seq, start = start_pos, end = stop_pos)

  # Translate the coding sequence into an amino acid sequence
  aa_seq <- translate_dna(coding_seq)

  # Write the amino acid sequence to a .txt file
  write(aa_seq, file = output_file)
}

####################################################################################################
# 2. All programmed codon variants (possible_mutations.csv)
####################################################################################################
# ------------------------------------------------------------------------------
# Script: Generate Programmed Codon Variants
# Description: Generates all possible programmed codon mutations for a given
#              wild-type sequence based on a specified mutagenesis strategy.
# Input:
#   - wt_seq_input: Wild-type sequence (string or path to FASTA file).
#   - start_stop_pos: Target sequence range format "start-stop".
#   - mutagenesis_type: Strategy ('nnk', 'nns', 'nnh', 'nnn', 'nnk_nns', 'nnk_nns_nnh', 'custom').
#   - custom_codon_library_path: Path to custom library. Automatically detects
#     if the file is a global list ("AAA, AAC...") or a position-wise CSV
#     (requires a "Position" header).
#   - output_file: Desired name/path for the output CSV.
# Output: A CSV file containing all possible mutated codons per position.
# ------------------------------------------------------------------------------

generate_possible_variants <- function(wt_seq_input, start_stop_pos, mutagenesis_type,
                                       custom_codon_library_path, output_file) {

  # Parse the start and stop positions from the input format "start-stop"
  positions <- unlist(strsplit(start_stop_pos, "-"))
  start_pos <- as.numeric(positions[1])
  stop_pos <- as.numeric(positions[2])

  # Load sequence from file or process as a direct string
  if (file.exists(wt_seq_input)) {
    wt_seq <- read_fasta(wt_seq_input)[1]
  } else {
    wt_seq <- toupper(wt_seq_input)
  }

  # Extract the target coding sequence
  coding_seq <- subseq_checked(wt_seq, start = start_pos, end = stop_pos)

  # Predefined codon dictionaries
  nnk_codons <- c('AAG', 'AAT', 'ATG', 'ATT', 'AGG', 'AGT', 'ACG', 'ACT', 'TAG', 'TAT', 'TTG', 'TTT', 'TGG', 'TGT', 'TCG', 'TCT', 'GAG', 'GAT', 'GTG', 'GTT', 'GGG', 'GGT', 'GCG', 'GCT', 'CAG', 'CAT', 'CTG', 'CTT', 'CGG', 'CGT', 'CCG', 'CCT')
  nns_codons <- c('AAG', 'AAC', 'ATG', 'ATC', 'AGG', 'AGC', 'ACG', 'ACC', 'TAG', 'TAC', 'TTG', 'TTC', 'TGG', 'TGC', 'TCG', 'TCC', 'GAG', 'GAC', 'GTG', 'GTC', 'GGG', 'GGC', 'GCG', 'GCC', 'CAG', 'CAC', 'CTG', 'CTC', 'CGG', 'CGC', 'CCG', 'CCC')
  nnh_codons <- c('AAA', 'AAC', 'AAT', 'ATA', 'ATC', 'ATT', 'AGA', 'AGC', 'AGT', 'ACA', 'ACC', 'ACT', 'TAA', 'TAC', 'TAT', 'TTA', 'TTC', 'TTT', 'TGA', 'TGC', 'TGT', 'TCA', 'TCC', 'TCT', 'GAA', 'GAC', 'GAT', 'GTA', 'GTC', 'GTT', 'GGA', 'GGC', 'GGT', 'GCA', 'GCC', 'GCT', 'CAA', 'CAC', 'CAT', 'CTA', 'CTC', 'CTT', 'CGA', 'CGC', 'CGT', 'CCA', 'CCC', 'CCT')
  nnn_codons <- c('AAA', 'AAC', 'AAG', 'AAT', 'ATA', 'ATC', 'ATG', 'ATT', 'AGA', 'AGC', 'AGG', 'AGT', 'ACA', 'ACC', 'ACG', 'ACT', 'TAA', 'TAC', 'TAG', 'TAT', 'TTA', 'TTC', 'TTG', 'TTT', 'TGA', 'TGC', 'TGG', 'TGT', 'TCA', 'TCC', 'TCG', 'TCT', 'GAA', 'GAC', 'GAG', 'GAT', 'GTA', 'GTC', 'GTG', 'GTT', 'GGA', 'GGC', 'GGG', 'GGT', 'GCA', 'GCC', 'GCG', 'GCT', 'CAA', 'CAC', 'CAG', 'CAT', 'CTA', 'CTC', 'CTG', 'CTT', 'CGA', 'CGC', 'CGG', 'CGT', 'CCA', 'CCC', 'CCG', 'CCT')

  # --------------------------------------------------------------------------
  # Custom Library Parsing with Auto-Detection
  # --------------------------------------------------------------------------
  is_position_wise <- FALSE
  position_lookup <- list()
  custom_codons <- NULL

  if (mutagenesis_type == "custom") {
    if (is.null(custom_codon_library_path) || !nzchar(custom_codon_library_path) ||
        !file.exists(custom_codon_library_path)) {
      stop("--mutagenesis_type custom requires a valid --custom_codon_library file.")
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
    } else if (mutagenesis_type == "nnh") {
            return(nnh_codons)
    } else if (mutagenesis_type == "nnn") {
            return(nnn_codons)
    } else if (mutagenesis_type == "nnk_nns") {
      if (substr(wt_codon, 3, 3) == "T") return(nns_codons) else return(nnk_codons)
    } else if (mutagenesis_type == "nnk_nns_nnh") {
            if (substr(wt_codon, 3, 3) == "T") {
                return(nns_codons)
            } else if (substr(wt_codon, 3, 3) == "G"){
                return(nnh_codons)
            } else {
                return(nnk_codons)
            }
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
      stop("Invalid mutagenesis_type. Choose from 'nnk', 'nns', 'nnh', 'nnn', 'nnk_nns', 'nnk_nns_nnh', or 'custom'.")
    }
  }

  # Iterate over each wild-type codon to assign programmed variants. Collect one data.frame per
  # position in a list and bind once at the end: growing `result` with rbind() inside the loop
  # reallocates the whole table on every variant, which is O(n^2) and becomes the dominant cost
  # for a long ORF with a permissive library (nnn gives 63 variants per position).
  per_position <- vector("list", length(wt_codons))

  for (i in seq_along(wt_codons)) {
    wt_codon <- wt_codons[i]
    codon_list <- get_codon_list(wt_codon, i)

    # Skip iteration if no custom codons were assigned to this specific position
    if (is.null(codon_list)) next

    # Filter out the wild-type codon from the mutation list
    possible_variants <- codon_list[codon_list != wt_codon]
    if (length(possible_variants) == 0) next

    per_position[[i]] <- data.frame(
      Codon_Number = i,
      wt_codon = wt_codon,
      Variant = possible_variants,
      stringsAsFactors = FALSE
    )
  }

  result <- do.call(rbind, per_position)
  if (is.null(result)) {
    result <- data.frame(Codon_Number = integer(0), wt_codon = character(0),
                         Variant = character(0), stringsAsFactors = FALSE)
  }

  write.csv(result, output_file, row.names = FALSE)
}

#####
# run functions
#####
aa_seq(
  wt_seq_input = "$wt_seq",
  "$pos_range",
  "aa_seq.txt"
  )

# Without '--custom_codon_library', Nextflow stages no file and the path is empty.
custom_lib_arg <- if ("$custom_codon_library" == "") {
  NULL
} else {
  "$custom_codon_library"
}

generate_possible_variants(
  wt_seq_input              = "$wt_seq",
  start_stop_pos            = "$pos_range",
  mutagenesis_type          = "$mutagenesis_type",
  custom_codon_library_path = custom_lib_arg,
  output_file               = "possible_mutations.csv"
)

####
# create versions.yml
####
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
