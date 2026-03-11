#!/usr/bin/env Rscript

suppressMessages(library(Biostrings))

generate_fitness_input <- function(wt_seq_path, gatk_file, pos_range, output_file_path) {
  # Parse the position range
  positions <- unlist(strsplit(pos_range, "-"))
  start_pos <- as.numeric(positions[1])
  stop_pos <- as.numeric(positions[2])

  # Load the wild-type sequence
  seq_data <- Biostrings::readDNAStringSet(filepath = wt_seq_path)
  wt_seq <- seq_data[[1]]  # Extract the sequence
  wt_seq <- subseq(wt_seq, start = start_pos, end = stop_pos)

  # Convert wt_seq to a character string
  wt_seq <- as.character(wt_seq)

  # Split the wild-type sequence into codons (groups of 3 bases)
  wt_codons <- substring(wt_seq, seq(1, nchar(wt_seq), 3), seq(3, nchar(wt_seq), 3))

  # Helper function to process GATK CSVs into count data
  process_gatk_file <- function(gatk_csv) {
    # Load the input GATK CSV file
    gatk_data <- read.csv(gatk_csv, stringsAsFactors = FALSE)

    # Initialize a data frame for results
    results <- data.frame(
      nt_seq = character(),
      count = numeric(),
      stringsAsFactors = FALSE
    )

    # Iterate over each row in the input data
    for (i in 1:nrow(gatk_data)) {
      # Extract the mutation info
      codon_mut <- gatk_data\$codon_mut[i]
      counts <- gatk_data\$counts[i]

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
      mutated_seq_string <- paste(mutated_codons, collapse = "")

      # Add the result to the data frame
      results <- rbind(results, data.frame(nt_seq = mutated_seq_string, count = counts))
    }

    return(results)
  }

  # Process the GATK file
  cat("Processing GATK file...\\n")
  processed_data <- process_gatk_file(gatk_file)

  # Write the processed data to a file without column names
  write.table(processed_data, file = output_file_path, sep = "\\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

#####
# run function
#####
generate_fitness_input(
  wt_seq_path = "$wt_seq",
  gatk_file = "$variantCounts_filtered_by_library",
  pos_range = "$pos_range",
  output_file_path = "${meta.id}_fitness_input.tsv"
)

#####
# create versions.yml
#####
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
