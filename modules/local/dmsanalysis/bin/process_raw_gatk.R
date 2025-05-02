# input: gatk variantcounts tsv file, output_path
# output: csv with column names. Creates additional counts_per_cov column. Fills pos_mut column for synonymous mutations. Sorted out variants that have mutations, but do not show up in the specifying columns -> was affecting roughly 30 low-count variants out of over 15000 in Taylor's data

library("dplyr")
process_raw_gatk <- function(gatk_file_path, output_csv_path) {

  # Set the column names
  colnames <- c("counts", "cov", "mean_length_variant_reads", "varying_bases",
                          "base_mut", "varying_codons", "codon_mut", "aa_mut", "pos_mut")

  # Read the GATK file into a data frame
  gatk_raw <- read.table(gatk_file_path, sep = "\t", header = FALSE, fill = TRUE, col.names = colnames)

  # Filter out rows where 'aa_mut' is empty or NA
  gatk_raw <- gatk_raw[!(gatk_raw$aa_mut == "" | is.na(gatk_raw$aa_mut)), ]

  # Handle synonymous mutations: where aa_mut starts with "S" and pos_mut is either NA or ""
  gatk_raw <- gatk_raw %>%
    rowwise() %>%
    mutate(
      pos_mut = ifelse(
        (is.na(pos_mut) | pos_mut == "") & grepl("^S:", aa_mut),
        # Construct the new 'pos_mut' entry for synonymous mutations
        paste0(
          sub("S:([A-Z])>[A-Z]", "\\1", aa_mut),  # Get the original amino acid from 'aa_mut'
          sub("^(\\d+):.*", "\\1", codon_mut),    # Get the position from 'codon_mut'
          sub("S:[A-Z]>([A-Z])", "\\1", aa_mut)   # Get the mutated amino acid from 'aa_mut'
        ),
        pos_mut  # Keep the existing 'pos_mut' if it's not NA or ""
      )
    ) %>%
    ungroup() %>%
    mutate(counts_per_cov = counts / cov)

  # Write the cleaned data frame to a CSV file
  write.csv(gatk_raw, file = output_csv_path, row.names = FALSE)
}




# Example usage (can be used for testing):
# process_raw_gatk("/path/to/gatk_file.txt", "/path/to/output_file.csv")
#process_raw_gatk("/Users/benjaminwehnert/CRG/DMS_QC/testing_data/output_premerged_vsearch.variantCounts", "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/raw_gatk.csv")


