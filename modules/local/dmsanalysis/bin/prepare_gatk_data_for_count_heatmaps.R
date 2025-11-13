# input: prefiltered GATK path (filtered for codon library), aa-seq file path, output path, threshold (for minimum counts to recognize variant)
# output: csv file serving as basis for counts_per_cov_heatmap function

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(reshape2))
suppressMessages(library(scales))

prepare_gatk_data_for_counts_heatmaps <- function(gatk_file_path, aa_seq_file_path, output_csv_path, threshold = 3) {
  # Load the raw GATK data
  raw_gatk <- read.table(gatk_file_path, sep = ",", header = TRUE)

  # Read the wild-type amino acid sequence from the text file
  wt_seq <- readLines(aa_seq_file_path)
  wt_seq <- unlist(strsplit(wt_seq, ""))  # Split the sequence into individual amino acids

  # Summarize counts-per-cov for each unique aa mutation in pos_mut
  aggregated_data <- raw_gatk %>%
    group_by(pos_mut) %>%
    summarize(total_counts_per_cov = sum(counts_per_cov, na.rm = TRUE),
              total_counts = sum(counts, na.rm = TRUE))  # Also sum the counts

  # Extract the wild-type position and mutations from 'pos_mut'
  aggregated_data <- aggregated_data %>%
    mutate(
      wt_aa = sub("(\\D)(\\d+)(\\D)", "\\1", pos_mut),  # Wild-type amino acid (e.g., S)
      position = as.numeric(sub("(\\D)(\\d+)(\\D)", "\\2", pos_mut)),  # Position (e.g., 3)
      mut_aa = sub("(\\D)(\\d+)(\\D)", "\\3", pos_mut)   # Mutant amino acid (e.g., R)
    )

  # Replace 'X' with '*', indicating the stop codon
  aggregated_data <- aggregated_data %>%
    mutate(mut_aa = ifelse(mut_aa == "X", "*", mut_aa))

  # Replace 'X' with '*' in the wild-type amino acid sequence as well
  wt_seq <- ifelse(wt_seq == "X", "*", wt_seq)

  # Define all 20 standard amino acids and the stop codon "*"
  all_amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                       "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "*")

  # Create a list of all positions in the wild-type sequence
  all_positions <- 1:length(wt_seq)

  # Create a complete grid of all possible combinations of positions and amino acids
  complete_data <- expand.grid(mut_aa = all_amino_acids, position = all_positions)

  # Merge the summarized data with the complete grid (filling missing entries with 0)
  heatmap_data <- complete_data %>%
    left_join(aggregated_data, by = c("mut_aa", "position")) %>%
    mutate(total_counts_per_cov = ifelse(is.na(total_counts_per_cov), 0, total_counts_per_cov),
           wt_aa = wt_seq[position])  # Assign the wild-type amino acid

  # Set variants with counts < threshold to NA
  heatmap_data <- heatmap_data %>%
    mutate(
      total_counts_per_cov = ifelse(total_counts < threshold, NA, total_counts_per_cov),
      total_counts = ifelse(total_counts < threshold, NA, total_counts)
    )

  # Fill pos_mut column
  heatmap_data <- heatmap_data %>%
    mutate(
      pos_mut = ifelse(is.na(pos_mut),
                       paste0(wt_aa, position, mut_aa),
                       pos_mut)
    )

  # Save the aggregated data to a CSV file
  write.csv(heatmap_data, file = output_csv_path, row.names = FALSE)
  print(paste("Aggregated data saved to:", output_csv_path))
}


# Aufruf der Datenaufbereitungsfunktion
# prepare_gatk_data_for_counts_heatmaps(
#   "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/gatk_filtered_by_codon_library.csv",
#   "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/aa_seq.txt",
#   "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/testing_outputs/prepared_gatk_data.csv",
#   threshold = 3
# )
