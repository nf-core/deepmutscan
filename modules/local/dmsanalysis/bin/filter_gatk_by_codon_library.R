# Input: pre_processed_raw_gatk_path, mutation library .csv path (former "possible_NNK_mutations.csv"), output_path
# Output: gatk table filtered for only single-codon mutations that are part of the library

# Load necessary library
library("dplyr")
library("stringr")

# filter_gatk_by_codon_library <- function(gatk_file_path, codon_library, output_file_path) {
#   # Load the GATK table from the provided file path
#   gatk_table <- read.csv(gatk_file_path)
#
#   # Predefined codon libraries
#   nnk_codons <- c('AAG', 'AAT', 'ATG', 'ATT', 'AGG', 'AGT', 'ACG', 'ACT',
#                   'TAG', 'TAT', 'TTG', 'TTT', 'TGG', 'TGT', 'TCG', 'TCT',
#                   'GAG', 'GAT', 'GTG', 'GTT', 'GGG', 'GGT', 'GCG', 'GCT',
#                   'CAG', 'CAT', 'CTG', 'CTT', 'CGG', 'CGT', 'CCG', 'CCT')
#
#   # Check if the codon_library is a predefined string or a custom vector
#   if (is.character(codon_library) && length(codon_library) == 1) {
#     if (codon_library == "NNK") {
#       codon_set <- nnk_codons
#     } else {
#       stop("Invalid predefined codon library specified.")
#     }
#   } else if (is.vector(codon_library)) {
#     codon_set <- codon_library
#   } else {
#     stop("Invalid codon library format. Must be a predefined string or a custom vector.")
#   }
#
#   # Filter for single-codon mutations that are part of the codon library and make sure to handle some mistaken formatting from gatk (there are cases where more than 3 Bases are mutated, but column varying_codons == 1)
#   filtered_gatk <- gatk_table %>%
#     filter(varying_codons == 1 & sub(".*>", "", codon_mut) %in% codon_set) %>%
#     rowwise() %>%
#     filter({
#       # Split base_mut into individual mutations
#       mutations <- unlist(strsplit(base_mut, ",\\s*"))  # Splits by comma and removes extra spaces
#       # Extract numeric positions from each mutation string
#       positions <- as.numeric(str_extract(mutations, "^[0-9]+"))
#
#       # Calculate the distance between the first and last position
#       distance <- max(positions, na.rm = TRUE) - min(positions, na.rm = TRUE)
#
#       # Keep rows where the distance is <= 2
#       distance <= 2
#     }) %>%
#     ungroup()
#
#
#   # Write the filtered GATK table to the output file path
#   write.csv(filtered_gatk, file = output_file_path, row.names = FALSE)
# }
#
# # Example usage with predefined library (NNK):
# # filtered_gatk <- filter_gatk_by_codon_library(raw_gatk, codon_library = "NNK")
#
# # Example usage with custom library:
# # custom_library <- c('AAG', 'TGT', 'GCT')
# # filtered_gatk <- filter_gatk_by_codon_library(raw_gatk, custom_library = custom_library)
#
# #filter_gatk_by_codon_library("/Users/benjaminwehnert/CRG/DMS_QC/testing_data/raw_gatk.csv", codon_library = "NNK", output_file_path = "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/gatk_filtered_by_codon_library.csv")





library("dplyr")
library("stringr")

filter_gatk_by_codon_library <- function(gatk_file_path, codon_library_path, output_file_path) {
  # Load the GATK table from the provided file path
  gatk_table <- read.csv(gatk_file_path)

  # Load the codon library from the provided .csv file
  codon_library <- read.csv(codon_library_path)

  # Ensure the codon library has the expected columns
  if (!all(c("Codon_Number", "wt_codon", "Variant") %in% colnames(codon_library))) {
    stop("Codon library must contain columns 'Codon_Number', 'wt_codon', and 'Variant'.")
  }

  # Filter the GATK table
  filtered_gatk <- gatk_table %>%
    filter(varying_codons == 1) %>% # Keep rows with single-codon mutations
    rowwise() %>%
    filter({
      # Extract the position and mutated codon
      codon_position <- as.numeric(sub(":.*", "", codon_mut)) # Extract position before ':'
      mutated_codon <- sub(".*>", "", codon_mut)             # Extract codon after '>'

      # Check if the position and codon are valid
      is_in_library <- any(
        codon_library$Codon_Number == codon_position &
          (codon_library$Variant == mutated_codon |  # Check Variant column
             codon_library$wt_codon == mutated_codon)  # Check wt_codon column
      )
      is_in_library
    }) %>%
    ungroup() %>%
    # Apply additional filtering based on mutation distances
    rowwise() %>%
    filter({
      # Split base_mut into individual mutations
      mutations <- unlist(strsplit(base_mut, ",\\s*"))  # Splits by comma and removes extra spaces
      # Extract numeric positions from each mutation string
      positions <- as.numeric(str_extract(mutations, "^[0-9]+"))

      # Calculate the distance between the first and last position
      distance <- max(positions, na.rm = TRUE) - min(positions, na.rm = TRUE)

      # Keep rows where the distance is <= 2
      distance <= 2
    }) %>%
    ungroup()

  # Write the filtered GATK table to the output file path
  write.csv(filtered_gatk, file = output_file_path, row.names = FALSE)
}

# example
#filter_gatk_by_codon_library("/Users/benjaminwehnert/CRG/DMS_QC/testing_data/raw_gatk.csv", codon_library = "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/possible_NNK_mutations.csv", output_file_path = "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/gatk_filtered_by_codon_library.csv")    ### this one's is correct for this data set
#filter_gatk_by_codon_library("/Users/benjaminwehnert/CRG/DMS_QC/testing_data/raw_gatk.csv", codon_library = "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/testing_outputs/possible_NNK_mutations_taylors_nnk_and_nns.csv", output_file_path = "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/testing_outputs/gatk_filtered_by_complete_codon_library.csv")


