# input: wt_seq_path, input-gatk-counts filtered by codon library, output-gatk-counts filtered by codon library, start-stop-codon
# output: merged counts file for input and output data that can be fed into DiMSum
# sums up all synonymous mutations, assigning them to the wildtype sequence in the first line


# suppressMessages(library(seqinr))
# suppressMessages(library(Biostrings))

# generate_dimsum_input <- function(wt_seq_path, gatk_input, gatk_output, pos_range) {
#   # Parse the position range
#   positions <- unlist(strsplit(pos_range, "-"))
#   start_pos <- as.numeric(positions[1])
#   stop_pos <- as.numeric(positions[2])
#
#   # Load the wild-type sequence
#   seq_data <- Biostrings::readDNAStringSet(filepath = wt_seq_path)
#   wt_seq <- seq_data[[1]]  # Extract the sequence
#   wt_seq <- subseq(wt_seq, start = start_pos, end = stop_pos)
#
#   # Convert wt_seq to a character string
#   wt_seq <- as.character(wt_seq)
#
#   # Split the wild-type sequence into codons (groups of 3 bases)
#   wt_codons <- substring(wt_seq, seq(1, nchar(wt_seq), 3), seq(3, nchar(wt_seq), 3))
#
#   # Helper function to process GATK CSVs into count data
#   process_gatk_file <- function(gatk_csv) {
#     # Load the input GATK CSV file
#     gatk_data <- read.csv(gatk_csv, stringsAsFactors = FALSE)
#
#     # Calculate the sum of all synonymous mutation counts
#     synonymous_counts <- sum(gatk_data$counts[grep("^S:", gatk_data$aa_mut)])
#
#     # Initialize a data frame with the wild-type sequence and its synonymous counts
#     results <- data.frame(
#       nt_seq = wt_seq,
#       count = synonymous_counts,
#       stringsAsFactors = FALSE
#     )
#
#     # Iterate over each row in the input data
#     for (i in 1:nrow(gatk_data)) {
#       # Extract the mutation info
#       codon_mut <- gatk_data$codon_mut[i]
#       counts <- gatk_data$counts[i]
#
#       # Create a mutable copy of the wild-type codons
#       mutated_codons <- wt_codons
#
#       # Apply the mutation
#       mutations <- strsplit(codon_mut, ", ")[[1]]
#       for (mutation in mutations) {
#         codon_position <- as.numeric(sub(":.*", "", mutation))
#         new_codon <- sub(".*>", "", mutation)
#         # Replace the codon at the specified position
#         mutated_codons[codon_position] <- new_codon
#       }
#
#       # Convert the mutated codons back to a sequence string
#       mutated_seq_string <- paste(mutated_codons, collapse = "")
#
#       # Add the result to the data frame
#       results <- rbind(results, data.frame(nt_seq = mutated_seq_string, count = counts))
#     }
#
#     return(results)
#   }
#
#   # Process the GATK input and output files
#   cat("Processing GATK input file...\n")
#   input_data <- process_gatk_file(gatk_input)
#   colnames(input_data)[2] <- "input1"  # Rename count column
#
#   cat("Processing GATK output file...\n")
#   output_data <- process_gatk_file(gatk_output)
#   colnames(output_data)[2] <- "output1"  # Rename count column
#
#   # Merge the input and output data
#   merged_data <- merge(input_data, output_data, by = "nt_seq", all = TRUE)
#   merged_data[is.na(merged_data)] <- 0  # Replace NA with 0 for missing counts
#
#   # Write the merged data to a file
#   output_file <- "merged_counts_for_dimsum.txt"
#   write.table(merged_data, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#   cat("Merged counts file created:", output_file, "\n")
# }
#

#
# generate_dimsum_input(
#   "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/MORtn5_reference.fa",
#   "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/gatk_filtered_by_codon_library.csv",
#   "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/gatk_filtered_by_codon_library_dummy.csv",
#   "23-1225"
# )


suppressMessages(library(Biostrings))

generate_dimsum_input <- function(wt_seq_path, gatk_file, pos_range, output_file_path) {
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
      codon_mut <- gatk_data$codon_mut[i]
      counts <- gatk_data$counts[i]

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
  cat("Processing GATK file...\n")
  processed_data <- process_gatk_file(gatk_file)

  # Write the processed data to a file without column names
  write.table(processed_data, file = output_file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# generate_dimsum_input(
#   "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/MORtn5_reference.fa",
#   "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/gatk_filtered_by_codon_library.csv",
#   "23-1225",
#   "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/testing_outputs/dimsum_input.tsv"
# )










