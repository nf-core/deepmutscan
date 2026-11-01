#!/usr/bin/env Rscript

####################################################################################################
# Base-R replacements for Biostrings (readDNAStringSet, subseq), with Biostrings' results
####################################################################################################

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

generate_fitness_input <- function(wt_seq_path, counts_file, pos_range, output_file_path) {
  # Parse the position range
  positions <- unlist(strsplit(pos_range, "-"))
  start_pos <- as.numeric(positions[1])
  stop_pos <- as.numeric(positions[2])

  # Load the wild-type sequence
  wt_seq <- read_fasta(wt_seq_path)[1]  # the first sequence
  wt_seq <- subseq_checked(wt_seq, start = start_pos, end = stop_pos)

  # Split the wild-type sequence into codons (groups of 3 bases)
  wt_codons <- substring(wt_seq, seq(1, nchar(wt_seq), 3), seq(3, nchar(wt_seq), 3))

  # Helper function to process variant-count CSVs into count data
  process_counts_file <- function(counts_csv) {
    # Load the input variant-count CSV file
    counts_data <- read.csv(counts_csv, stringsAsFactors = FALSE)

    # Pre-allocate one slot per row. Growing `results` with rbind() inside the loop reallocated
    # the whole frame on every variant (O(n^2)); for a long ORF that dominated the runtime.
    n_rows <- nrow(counts_data)
    seq_strings <- character(n_rows)
    seq_counts  <- numeric(n_rows)

    if (n_rows == 0) {
      return(data.frame(nt_seq = character(0), count = numeric(0), stringsAsFactors = FALSE))
    }

    # Iterate over each row in the input data
    for (i in 1:n_rows) {
      # Extract the mutation info
      codon_mut <- counts_data\$codon_mut[i]

      if("counts_corrected" %in% colnames(counts_data)){
        counts <- counts_data\$counts_corrected[i]
      }else{
        counts <- counts_data\$counts[i]
      }

      # Create a mutable copy of the wild-type codons
      mutated_codons <- wt_codons

      # Apply the mutation
      mutations <- strsplit(codon_mut, ", ")[[1]]
      for (mutation in mutations) {
        codon_position <- as.numeric(sub(":.*", "", mutation))
        new_codon <- sub(".*>", "", mutation)
        # Replace the codon at the specified position
        mutated_codons[codon_position] <- new_codon
      }

      # Convert the mutated codons back to a sequence string
      seq_strings[i] <- paste(mutated_codons, collapse = "")
      seq_counts[i]  <- counts
    }

    return(data.frame(nt_seq = seq_strings, count = seq_counts, stringsAsFactors = FALSE))
  }

  # Process the variant-count file
  cat("Processing variant-count file...\\n")
  processed_data <- process_counts_file(counts_file)

  # Write the processed data to a file without column names
  write.table(processed_data, file = output_file_path, sep = "\\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

#####
# run function
#####
generate_fitness_input(
  wt_seq_path = "$wt_seq",
  counts_file = "$variantCounts_filtered_by_library",
  pos_range = "$pos_range",
  output_file_path = "${meta.id}_fitness_input.tsv"
)

#####
# create versions.yml
#####
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
