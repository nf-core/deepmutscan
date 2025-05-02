# input: NNK_codon_library_filtered_gatk.csv-path, prefiltered_gatk.csv-path (containing only NNK mutations), output_folder-path
# output: completed gatk_file with all possible variants (even if not measured in sequencing) -> NA in counts and counts_per_cov to 0.0000001 to deal with log-scale in following calculations

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

complete_prefiltered_gatk <- function(possible_nnk_path, prefiltered_gatk_path, output_file_path) {

  # Load the possible NNK mutations CSV
  possible_nnk <- read.csv(possible_nnk_path)

  # Load the prefiltered GATK CSV
  prefiltered_gatk <- read.csv(prefiltered_gatk_path)

  # Create codon_mut column in possible_NNK_mutations in the format 'Codon_Number:wt_codon>Variant'
  possible_nnk <- possible_nnk %>%
    mutate(codon_mut = paste0(Codon_Number, ":", wt_codon, ">", Variant))

  # Merge both dataframes based on the codon_mut column (full join to include all)
  merged_data <- full_join(prefiltered_gatk, possible_nnk, by = "codon_mut")

  # Fill missing values in counts_per_cov and counts with 0.0000001
  merged_data <- merged_data %>%
    mutate(counts_per_cov = ifelse(is.na(counts_per_cov), 0.0000001, counts_per_cov),
           counts = ifelse(is.na(counts), 0.000001, counts))

  # Calculate Hamming distance (varying_bases) and mutation details (aa_mut, pos_mut)
  merged_data <- merged_data %>%
    rowwise() %>%
    mutate(varying_bases = hamming_distance(wt_codon, Variant),
           mutation_info = list(mutation_details(wt_codon, Variant, Codon_Number))) %>%
    mutate(aa_mut = mutation_info$aa_mut,  # Extract aa_mut
           pos_mut = mutation_info$pos_mut) %>%  # Extract pos_mut
    ungroup() %>%
    select(-mutation_info)  # Remove the temporary list column

  # Save the merged data to a new CSV file
  write.csv(merged_data, file = output_file_path, row.names = FALSE)
}

# Example call
#complete_prefiltered_gatk("/Users/benjaminwehnert/CRG/DMS_QC/testing_data/possible_NNK_mutations.csv", "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/gatk_filtered_by_codon_library.csv", "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/testing_outputs")

