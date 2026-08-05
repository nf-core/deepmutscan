#!/usr/bin/env Rscript

# four parts condensed into one R script

#########
# process raw counts
#########

# input: counts variantcounts tsv file, output_path
# output: csv with column names. Creates additional counts_per_cov column. Fills pos_mut column for synonymous mutations. Sorted out variants that have mutations, but do not show up in the specifying columns -> was affecting roughly 30 low-count variants out of over 15000 in Taylor's data

library("dplyr")
process_raw_counts <- function(counts_file_path, output_csv_path) {

  # Set the column names
  colnames <- c("counts", "cov", "mean_length_variant_reads", "varying_bases",
                          "base_mut", "varying_codons", "codon_mut", "aa_mut", "pos_mut")

  # Read the variant-count file into a data frame
  counts_raw <- read.table(counts_file_path, sep = "\\t", header = FALSE, fill = TRUE, col.names = colnames)

  # Filter out rows where 'aa_mut' is empty or NA
  counts_raw <- counts_raw[!(counts_raw\$aa_mut == "" | is.na(counts_raw\$aa_mut)), ]

  # Handle synonymous mutations: where aa_mut starts with "S" and pos_mut is either NA or ""
  counts_raw <- counts_raw %>%
    rowwise() %>%
    mutate(
      pos_mut = ifelse(
        (is.na(pos_mut) | pos_mut == "") & grepl("^S:", aa_mut),
        # Construct the new 'pos_mut' entry for synonymous mutations
        paste0(
          sub("S:([A-Z])>[A-Z]", "\\\\1", aa_mut),  # Get the original amino acid from 'aa_mut'
          sub("^(\\\\d+):.*", "\\\\1", codon_mut),    # Get the position from 'codon_mut'
          sub("S:[A-Z]>([A-Z])", "\\\\1", aa_mut)   # Get the mutated amino acid from 'aa_mut'
        ),
        pos_mut  # Keep the existing 'pos_mut' if it's not NA or ""
      )
    ) %>%
    ungroup() %>%
    mutate(counts_per_cov = counts / cov)

  # Write the cleaned data frame to a CSV file
  write.csv(counts_raw, file = output_csv_path, row.names = FALSE)
}

#####
# run function
#####
process_raw_counts(
  counts_file_path = "$variantCounts",
  output_csv_path = "annotated_variantCounts.csv"
)








#########
# filter counts by codon library
#########

# Input: pre_processed_raw_counts_path, mutation library .csv path (former "possible_NNK_mutations.csv"), output_path
# Output: counts table filtered for only single-codon mutations that are part of the library

library("dplyr")
library("stringr")

filter_counts_by_codon_library <- function(counts_file_path, codon_library_path, output_file_path) {
  # Load the variant-count table from the provided file path
  counts_table <- read.csv(counts_file_path)

  # Load the codon library from the provided .csv file
  codon_library <- read.csv(codon_library_path)

  # Ensure the codon library has the expected columns
  if (!all(c("Codon_Number", "wt_codon", "Variant") %in% colnames(codon_library))) {
    stop("Codon library must contain columns 'Codon_Number', 'wt_codon', and 'Variant'.")
  }

  # Filter the variant-count table
  filtered_counts <- counts_table %>%
    filter(varying_codons == 1) %>% # Keep rows with single-codon mutations
    rowwise() %>%
    filter({
      # Extract the position and mutated codon
      codon_position <- as.numeric(sub(":.*", "", codon_mut)) # Extract position before ':'
      mutated_codon <- sub(".*>", "", codon_mut)             # Extract codon after '>'

      # Check if the position and codon are valid
      is_in_library <- any(
        codon_library\$Codon_Number == codon_position &
          (codon_library\$Variant == mutated_codon |  # Check Variant column
             codon_library\$wt_codon == mutated_codon)  # Check wt_codon column
      )
      is_in_library
    }) %>%
    ungroup() %>%
    # Apply additional filtering based on mutation distances
    rowwise() %>%
    filter({
      # Split base_mut into individual mutations
      mutations <- unlist(strsplit(base_mut, ",\\\\s*"))  # Splits by comma and removes extra spaces
      # Extract numeric positions from each mutation string
      positions <- as.numeric(str_extract(mutations, "^[0-9]+"))

      # Calculate the distance between the first and last position
      distance <- max(positions, na.rm = TRUE) - min(positions, na.rm = TRUE)

      # Keep rows where the distance is <= 2
      distance <= 2
    }) %>%
    ungroup()

  # Write the filtered variant-count table to the output file path
  write.csv(filtered_counts, file = output_file_path, row.names = FALSE)
}

#####
# run function
#####
filter_counts_by_codon_library(
  counts_file_path = "annotated_variantCounts.csv",
  codon_library_path = "$possible_mutations",
  output_file_path = "variantCounts_filtered_by_library.csv"
)




#########
# complete filtered counts
#########

# input: NNK_codon_library_filtered_counts.csv-path, prefiltered_counts.csv-path (containing only NNK mutations), output_folder-path
# output: completed counts_file with all possible variants (even if not measured in sequencing) -> NA in counts and counts_per_cov to 0.0000001 to deal with log-scale in following calculations

library(dplyr)
library(Biostrings)  # Required for codon-to-amino-acid translation

# Function to calculate Hamming distance (varying_bases)
hamming_distance <- function(wt_codon, variant_codon) {
  sum(strsplit(wt_codon, "")[[1]] != strsplit(variant_codon, "")[[1]])
}

# Function to get amino acid from codon
get_amino_acid <- function(codon) {
  codon_table <- GENETIC_CODE
  aa <- codon_table[[toupper(codon)]]
  if (is.null(aa)) {
    return(NA)  # Handle cases where codon is not valid
  }
  return(aa)
}

# Function to calculate mutation type (aa_mut) and pos_mut
mutation_details <- function(wt_codon, variant_codon, codon_number) {
  wt_aa <- get_amino_acid(wt_codon)
  variant_aa <- get_amino_acid(variant_codon)

  # If amino acids are different, it's a missense mutation; otherwise, synonymous
  if (wt_aa != variant_aa) {
    mutation_type <- "M"  # Missense mutation
  } else {
    mutation_type <- "S"  # Synonymous mutation
  }

  # aa_mut: Type of mutation and amino acid changes (e.g., M:D>S)
  aa_mut <- paste0(mutation_type, ":", wt_aa, ">", variant_aa)

  # pos_mut: Wild-type AA, codon position, mutated AA (e.g., D2Q)
  pos_mut <- paste0(wt_aa, codon_number, variant_aa)

  return(list(aa_mut = aa_mut, pos_mut = pos_mut))
}

complete_prefiltered_counts <- function(possible_nnk_path, prefiltered_counts_path, output_file_path) {

  # Load the possible NNK mutations CSV
  possible_nnk <- read.csv(possible_nnk_path)

  # Load the prefiltered variant-count CSV
  prefiltered_counts <- read.csv(prefiltered_counts_path)

  # Create codon_mut column in possible_NNK_mutations in the format 'Codon_Number:wt_codon>Variant'
  possible_nnk <- possible_nnk %>%
    mutate(codon_mut = paste0(Codon_Number, ":", wt_codon, ">", Variant))

  # Merge both dataframes based on the codon_mut column (full join to include all)
  prefiltered_counts\$codon_mut <- toupper(prefiltered_counts\$codon_mut)
  merged_data <- full_join(prefiltered_counts, possible_nnk, by = "codon_mut")

  # Fill missing values in counts_per_cov and counts with 0.0000001
  merged_data <- merged_data %>%
    mutate(counts_per_cov = ifelse(is.na(counts_per_cov), 0.0000001, counts_per_cov),
           counts = ifelse(is.na(counts), 0.000001, counts))

  # Calculate Hamming distance (varying_bases) and mutation details (aa_mut, pos_mut)
  merged_data <- merged_data %>%
    rowwise() %>%
    mutate(varying_bases = hamming_distance(wt_codon, Variant),
           mutation_info = list(mutation_details(wt_codon, Variant, Codon_Number))) %>%
    mutate(aa_mut = mutation_info\$aa_mut,  # Extract aa_mut
           pos_mut = mutation_info\$pos_mut) %>%  # Extract pos_mut
    ungroup() %>%
    select(-mutation_info)  # Remove the temporary list column

  # Save the merged data to a new CSV file
  write.csv(merged_data, file = output_file_path, row.names = FALSE)
}

#####
# run function
#####
complete_prefiltered_counts(
  possible_nnk_path = "$possible_mutations",
  prefiltered_counts_path = "variantCounts_filtered_by_library.csv",
  output_file_path = "library_completed_variantCounts.csv"
)









#########
# prepare counts data for counts heatmap
#########

# input: prefiltered variant-count path (filtered for codon library), aa-seq file path, output path, threshold (for minimum counts to recognize variant)
# output: csv file serving as basis for counts_per_cov_heatmap function

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(reshape2))
suppressMessages(library(scales))

prepare_counts_data_for_counts_heatmaps <- function(counts_file_path, aa_seq_file_path, output_csv_path, threshold = 3) {
  # Load the raw variant-count data
  raw_counts <- read.table(counts_file_path, sep = ",", header = TRUE)

  # Read the wild-type amino acid sequence from the text file
  wt_seq <- readLines(aa_seq_file_path)
  wt_seq <- unlist(strsplit(wt_seq, ""))  # Split the sequence into individual amino acids

  # Summarize counts-per-cov for each unique aa mutation in pos_mut
  aggregated_data <- raw_counts %>%
    group_by(pos_mut) %>%
    summarize(total_counts_per_cov = sum(counts_per_cov, na.rm = TRUE),
              total_counts = sum(counts, na.rm = TRUE))  # Also sum the counts

  # Extract the wild-type position and mutations from 'pos_mut'
  aggregated_data <- aggregated_data %>%
    mutate(
      wt_aa = sub("(\\\\D)(\\\\d+)(\\\\D)", "\\\\1", pos_mut),  # Wild-type amino acid (e.g., S)
      position = as.numeric(sub("(\\\\D)(\\\\d+)(\\\\D)", "\\\\2", pos_mut)),  # Position (e.g., 3)
      mut_aa = sub("(\\\\D)(\\\\d+)(\\\\D)", "\\\\3", pos_mut)   # Mutant amino acid (e.g., R)
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


#####
# run function
#####
prepare_counts_data_for_counts_heatmaps(
  counts_file_path = "variantCounts_filtered_by_library.csv",
  aa_seq_file_path = "$aa_seq",
  output_csv_path = "variantCounts_for_heatmaps.csv",
  threshold = $min_counts
)







####
# create versions.yml
####
r_version <- strsplit(version[['version.string']], ' ')[[1]][3]
dplyr_version <- as.character(packageVersion("dplyr"))
Biostrings_version <- as.character(packageVersion("Biostrings"))
stringr_version <- as.character(packageVersion("stringr"))
ggplot2_version <- as.character(packageVersion("ggplot2"))
tidyr_version <- as.character(packageVersion("tidyr"))
reshape2_version <- as.character(packageVersion("reshape2"))
scales_version <- as.character(packageVersion("scales"))


if (is.null(r_version)) r_version <- "unknown"
if (length(dplyr_version) == 0) dplyr_version <- "unknown"
if (length(Biostrings_version) == 0) Biostrings_version <- "unknown"
if (length(stringr_version) == 0) stringr_version <- "unknown"
if (length(ggplot2_version) == 0) ggplot2_version <- "unknown"
if (length(tidyr_version) == 0) tidyr_version <- "unknown"
if (length(reshape2_version) == 0) reshape2_version <- "unknown"
if (length(scales_version) == 0) scales_version <- "unknown"


f <- file("versions.yml", "w")
writeLines(
  c(
    '"${task.process}":',
    paste('    r-base:', r_version),
    paste('    r-dplyr:', dplyr_version),
    paste('    r-Biostrings:', Biostrings_version),
    paste('    r-stringr:', stringr_version),
    paste('    r-ggplot2:', ggplot2_version),
    paste('    r-tidyr:', tidyr_version),
    paste('    r-reshape2:', reshape2_version),
    paste('    r-scales:', scales_version)
  ),
  f
)
close(f)
