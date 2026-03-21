#!/usr/bin/env Rscript

## WT-based sequencing error correction in nf-core/deepmutscan
## 21.03.2026


## --- Helper functiosn ---

# function to run the WT sequencing error correction (nucleotide level)
seq_error_correct_by_WT_nt <- function(wt_seq_count_path, input_count_path, output_file_path){
  
  ## load data (nucleotide-level counts from GATK)
  WT.counts <- read.csv(wt_seq_count_path)
  input.counts <- read.csv(input_count_path)
  
  ## append key columns
  input.counts$counts_corrected <- input.counts$counts
  input.counts$counts_per_cov_corrected <- input.counts$counts_per_cov
  
  # Process the GATK file
  for (i in 1:nrow(WT.counts)){
    
    ## only look at variants observed in both the input and WT sequencing
    tmp.id <- match(WT.counts[i,"base_mut"], input.counts[,"base_mut"])
    if(is.na(tmp.id) == T){
      next
    }
    
    ## subtract the observed per-base coverage in the WT sequencing
    input.counts[tmp.id,"counts_per_cov_corrected"] <- input.counts[tmp.id,"counts_per_cov"] - WT.counts[i,"counts_per_cov"]
    
    ## based on this, also adjust the total_counts (cross-multiplication)
    ## total_counts_corrected ~ total_counts * c(total_counts_per_cov_corrected / total_counts_per_cov)
    input.counts[tmp.id,"counts_corrected"] <- input.counts[tmp.id,"counts"] * c(input.counts[tmp.id,"counts_per_cov_corrected"] / input.counts[tmp.id,"counts_per_cov"])
    
    ## round to nearest integer, do not allow for negative counts
    input.counts[tmp.id,"counts_corrected"] <- round(input.counts[tmp.id,"counts_corrected"])
    if(input.counts[tmp.id,"counts_corrected"] < 0){
      input.counts[tmp.id,"counts_corrected"] <- 0
      input.counts[tmp.id,"counts_per_cov_corrected"] <- 0
    }
    
  }

  # Write the processed data
  write.csv(input.counts, file = output_file_path)
  
}

# to re-make the AA level input table after WT sequencing error correction
seq_error_correct_counts_for_heatmaps <- function(gatk_file_path, aa_seq_file_path, output_csv_path, threshold = 3) {
  
  # Load the raw GATK data
  raw_gatk <- read.table(gatk_file_path, sep = ",", header = TRUE, stringsAsFactors = FALSE)
  
  # Read the wild-type amino acid sequence from the text file
  wt_seq <- readLines(aa_seq_file_path)
  wt_seq <- unlist(strsplit(wt_seq, ""))
  
  # Replace 'X' with '*', indicating the stop codon
  wt_seq[wt_seq == "X"] <- "*"
  
  # Summarize counts for each unique pos_mut
  aggregated_counts_per_cov <- aggregate(
    counts_per_cov_corrected ~ pos_mut,
    data = raw_gatk,
    FUN = function(x) sum(x, na.rm = TRUE)
  )
  
  aggregated_counts <- aggregate(
    counts_corrected ~ pos_mut,
    data = raw_gatk,
    FUN = function(x) sum(x, na.rm = TRUE)
  )
  
  # Merge the two aggregated tables
  aggregated_data <- merge(
    aggregated_counts_per_cov,
    aggregated_counts,
    by = "pos_mut",
    all = TRUE
  )
  
  # Rename columns to match original output
  names(aggregated_data)[names(aggregated_data) == "counts_per_cov_corrected"] <- "total_counts_per_cov_corrected"
  names(aggregated_data)[names(aggregated_data) == "counts_corrected"] <- "total_counts_corrected"
  
  # Extract wt_aa, position, and mut_aa from pos_mut
  aggregated_data$wt_aa <- sub("(\\D)(\\d+)(\\D)", "\\1", aggregated_data$pos_mut)
  aggregated_data$position <- as.numeric(sub("(\\D)(\\d+)(\\D)", "\\2", aggregated_data$pos_mut))
  aggregated_data$mut_aa <- sub("(\\D)(\\d+)(\\D)", "\\3", aggregated_data$pos_mut)
  
  # Replace 'X' with '*'
  aggregated_data$mut_aa[aggregated_data$mut_aa == "X"] <- "*"
  
  # Define all 20 standard amino acids and stop codon
  all_amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                       "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "*")
  
  # Create all positions
  all_positions <- seq_along(wt_seq)
  
  # Create complete grid of all possible position/mutation combinations
  complete_data <- expand.grid(
    mut_aa = all_amino_acids,
    position = all_positions,
    stringsAsFactors = FALSE
  )
  
  # Merge aggregated data into complete grid
  heatmap_data <- merge(
    complete_data,
    aggregated_data,
    by = c("mut_aa", "position"),
    all.x = TRUE,
    sort = FALSE
  )
  
  # Fill missing values
  heatmap_data$total_counts_per_cov_corrected[is.na(heatmap_data$total_counts_per_cov_corrected)] <- 0
  heatmap_data$wt_aa <- wt_seq[heatmap_data$position]
  
  # Apply threshold
  low_count <- is.na(heatmap_data$total_counts_corrected) | heatmap_data$total_counts_corrected < threshold
  heatmap_data$total_counts_per_cov_corrected[low_count] <- NA
  heatmap_data$total_counts_corrected[low_count] <- NA
  
  # Fill missing pos_mut values
  missing_pos_mut <- is.na(heatmap_data$pos_mut)
  heatmap_data$pos_mut[missing_pos_mut] <- paste0(
    heatmap_data$wt_aa[missing_pos_mut],
    heatmap_data$position[missing_pos_mut],
    heatmap_data$mut_aa[missing_pos_mut]
  )
  
  # Re-order rows
  out.order <- paste0(rep(1:max(heatmap_data$position), each = 21),
                      rep(c("A", "C", "D", "E", "F", "G", "H",
                            "I", "K", "L","M", "N", "P", "Q",
                            "R", "S", "T", "V", "W", "Y", "*"),
                          max(heatmap_data$position)))
  heatmap_data <- heatmap_data[match(out.order, paste0(heatmap_data$position, heatmap_data$mut_aa)),]
  rownames(heatmap_data) <- 1:nrow(heatmap_data)
  
  # Save output
  write.csv(heatmap_data, file = output_csv_path, row.names = FALSE)

}

## --- Main functions ---
seq_error_correct_by_WT_nt(wt_seq_count_path = "intermediate_files/processed_gatk_files/wildtype_pe/variantCounts_filtered_by_library.csv",
                           input_count_path = "intermediate_files/processed_gatk_files/input_1_pe/variantCounts_filtered_by_library.csv",
                           output_file_path = "intermediate_files/processed_gatk_files/input_1_pe/variantCounts_filtered_by_library_err_corrected.csv")

seq_error_correct_counts_for_heatmaps(gatk_file_path = "intermediate_files/processed_gatk_files/input_1_pe/variantCounts_filtered_by_library_err_corrected.csv",
                                      aa_seq_file_path = "intermediate_files/aa_seq.txt",
                                      output_csv_path = "intermediate_files/processed_gatk_files/input_1_pe/variantCounts_for_heatmaps_err_corrected.csv")

####
# create versions.yml
####
r_version <- strsplit(version[['version.string']], ' ')[[1]][3]
if (is.null(r_version)) r_version <- "unknown"
f <- file("versions.yml", "w")
writeLines(
  c(
    '"\${task.process}":',
    paste('    r-base:', r_version),
  ),
  f
)
close(f)
