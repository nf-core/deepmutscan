suppressMessages(library(dplyr))

prepare_gatk_data_for_fitness_heatmap <- function(csv_file_path, aa_seq_path, dimsum_fitness_path, output_csv_path) {
  # Load the CSV file
  gatk_data <- read.csv(csv_file_path)

  # Read the wild-type amino acid sequence from the text file
  aa_seq <- readLines(aa_seq_path)
  aa_seq <- unlist(strsplit(aa_seq, ""))  # Split the sequence into individual amino acids

  # Add the mutated sequence column
  gatk_data <- gatk_data %>%
    mutate(
      mutated_sequence = sapply(1:nrow(gatk_data), function(i) {
        # Create a copy of the wild-type sequence
        mutated_seq <- aa_seq
        # Extract mutation information
        position <- as.numeric(position[i])
        mut_aa <- mut_aa[i]
        # Apply the mutation
        if (!is.na(mut_aa) && position > 0 && position <= length(mutated_seq)) {
          mutated_seq[position] <- mut_aa
        }
        # Return the mutated sequence as a string
        paste(mutated_seq, collapse = "")
      })
    )

  # # Save the updated data to a new CSV file
  # write.csv(gatk_data, file = output_csv_path, row.names = FALSE)
  # print(paste("Updated data with mutated sequences saved to:", output_csv_path))


  load(dimsum_fitness_path)
  dimsum_fitness <- rbind(all_variants, synonymous)
  dimsum_fitness <- dimsum_fitness[-which(dimsum_fitness$WT == T)[-1],] # remove duplicate
  rm(doubles,singles,all_variants,synonymous,wildtype)

  # Perform a left join to retain all rows from gatk_data
  merged_data <- gatk_data %>%
    left_join(dimsum_fitness, by = c("mutated_sequence" = "aa_seq"))

  # Handle wild-type rows specifically
  merged_data <- merged_data %>%
    mutate(
      # If WT from dimsum_fitness exists, use it; otherwise keep as is
      WT = ifelse(is.na(WT) & mutated_sequence == paste(aa_seq, collapse = ""), TRUE, WT)
    )

  # Save the merged data to a CSV file
  write.csv(merged_data, file = output_csv_path, row.names = FALSE)
  print(paste("Merged data saved to:", output_csv_path))
}




# add_mutated_sequence("/Users/benjaminwehnert/CRG/DMS_QC/testing_data/testing_outputs/prepared_gatk_data.csv",
#                      "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/aa_seq.txt",
#
#                      "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/testing_outputs/prepared_gatk_for_fitness.csv")






suppressMessages(library(dplyr))

add_mutated_sequence <- function(csv_file_path, aa_seq_path, dimsum_fitness_path, output_csv_path) {

  # Load the CSV file
  gatk_data <- read.csv(csv_file_path)

  # Read the wild-type amino acid sequence from the text file
  aa_seq <- readLines(aa_seq_path)
  aa_seq <- unlist(strsplit(aa_seq, ""))  # Split the sequence into individual amino acids

  # Add the mutated sequence column
  gatk_data <- gatk_data %>%
    mutate(
      mutated_sequence = sapply(1:nrow(gatk_data), function(i) {
        # Create a copy of the wild-type sequence
        mutated_seq <- aa_seq
        # Extract mutation information
        position <- as.numeric(gatk_data$position[i])
        mut_aa <- gatk_data$mut_aa[i]
        # Apply the mutation
        if (!is.na(mut_aa) && position > 0 && position <= length(mutated_seq)) {
          mutated_seq[position] <- mut_aa
        }
        # Return the mutated sequence as a string
        paste(mutated_seq, collapse = "")
      })
    )

  # Load dimsum_fitness data
  load(dimsum_fitness_path)
  dimsum_fitness <- rbind(all_variants, synonymous)

  # Remove unnecessary columns
  dimsum_fitness <- dimsum_fitness %>%
    select(-nt_seq, -Nham_nt, -Nmut_codons, -indel, -STOP, -STOP_readthrough)

  # Ensure one-to-one mapping for the merge
  dimsum_fitness <- dimsum_fitness %>%
    distinct(aa_seq, .keep_all = TRUE)  # Keep only one row per unique `aa_seq`

  # Perform a left join to merge gatk_data with dimsum_fitness
  merged_data <- gatk_data %>%
    left_join(dimsum_fitness, by = c("mutated_sequence" = "aa_seq"))

  # Save the merged data to a CSV file
  write.csv(merged_data, file = output_csv_path, row.names = FALSE)
  print(paste("Merged data saved to:", output_csv_path))
}







# test data
prepare_gatk_data_for_fitness_heatmap("/Users/benjaminwehnert/CRG/DMS_QC/howard/bin4/intermediate_files/processed_gatk_files/variantCounts_for_heatmaps.csv",
                     "/Users/benjaminwehnert/CRG/DMS_QC/howard/bin4/intermediate_files/aa_seq.txt",
                     "/Users/benjaminwehnert/CRG/DMS_QC/howard/dimsum/howards_data/howards_data_fitness_replicates.RData",
                     "/Users/benjaminwehnert/CRG/DMS_QC/howard/prepared_gatk_for_fitness.csv")

#test on gid1a data
prepare_gatk_data_for_fitness_heatmap("/Users/benjaminwehnert/Downloads/variantCounts_for_heatmaps.csv",
                                      "/Users/benjaminwehnert/Downloads/aa_seq.txt",
                                      "/Users/benjaminwehnert/Downloads/2025-01-20_GID1A_DiMSum_1_fitness_replicates.RData",
                                      "/Users/benjaminwehnert/Downloads/dimsum_fitness_gid1a/prepared_gatk_for_fitness.csv")
