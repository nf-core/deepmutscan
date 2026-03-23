#!/usr/bin/env Rscript

## false-doubles based sequencing error correction in nf-core/deepmutscan
## 23.03.2026


## --- Helper functions ---

# function to run the false-doubles based sequencing error correction (nucleotide level)
seq_error_correct_by_false_doubles <- function(input_count_path_raw, input_count_path_processed, output_file_path){
  
  ## load data (nucleotide-level counts from GATK)
  
  # Set the column names
  colnames <- c("counts", "cov", "mean_length_variant_reads", "varying_bases",
                "base_mut", "varying_codons", "codon_mut", "aa_mut", "pos_mut")
  input.counts.raw <- read.table(input_count_path_raw, sep = "\t", header = FALSE, fill = TRUE, col.names = colnames)
  input.counts.processed <- read.csv(input_count_path_processed)
  
  ## append key columns
  input.counts.processed$counts_corrected <- input.counts.processed$counts
  input.counts.processed$counts_per_cov_corrected <- input.counts.processed$counts_per_cov
  
  # Process the GATK file, only look at single nucleotide variants
  # Need to parallelise
  cat("Sequencing error correction of GATK counts...\n")
  for (i in grep("[,]", input.counts.processed[,"base_mut"], invert = T)){
    
    ## only look at single nucleotide variants
    
    ## locate this variant across all multi-codon variants in the original count matrix
    tmp.false.doubles <- input.counts.raw[grep(input.counts.processed[i,"codon_mut"], input.counts.raw[,"codon_mut"]),]
    
    ## subset multi-codon variants
    tmp.false.doubles <- tmp.false.doubles[which(tmp.false.doubles[,"varying_bases"] >= 3 & tmp.false.doubles[,"varying_codons"] < 3),]
    if(nrow(tmp.false.doubles) == 0){
      next
    }
    
    ## subtract the (probable) false fraction from actual single counts
    # boxplot(tmp.false.doubles[,"counts"] / tmp.false.doubles[,"cov"], log = "y", pch = 16)
    # median(tmp.false.doubles[,"counts"] / tmp.false.doubles[,"cov"])
    tmp.false.doubles.median.counts.per.cov <- median(tmp.false.doubles[,"counts"] / tmp.false.doubles[,"cov"])
    input.counts.processed[i,"counts_per_cov_corrected"] <- input.counts.processed[i,"counts_per_cov"] - sqrt(tmp.false.doubles.median.counts.per.cov) ## square root
    
    ## based on this, also adjust the total_counts (cross-multiplication)
    ## total_counts_corrected ~ total_counts * c(total_counts_per_cov_corrected / total_counts_per_cov)
    input.counts.processed[i,"counts_corrected"] <- input.counts.processed[i,"counts"] * c(input.counts.processed[i,"counts_per_cov_corrected"] / input.counts.processed[i,"counts_per_cov"])
    
    ## round to nearest integer, do not allow for negative counts
    input.counts.processed[i,"counts_corrected"] <- round(input.counts.processed[i,"counts_corrected"])
    if(input.counts.processed[i,"counts_corrected"] < 0){
      input.counts.processed[i,"counts_corrected"] <- 0
      input.counts.processed[i,"counts_per_cov_corrected"] <- 0
    }
    
  }
  
  # Write the processed data
  write.csv(input.counts.processed, file = output_file_path, row.names = F)
  
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
  print(paste("Aggregated data saved to:", output_csv_path))

}

## --- Main functions ---
seq_error_correct_by_false_doubles(input_count_path_raw = "intermediate_files/gatk/input_1_pe/gatk_output.variantCounts",
                                   input_count_path_processed = "intermediate_files/processed_gatk_files/input_1_pe/variantCounts_filtered_by_library.csv",
                                   output_file_path = "intermediate_files/processed_gatk_files/input_1_pe/variantCounts_filtered_by_library_err_corrected_false_doubles.csv")

seq_error_correct_counts_for_heatmaps(gatk_file_path = "intermediate_files/processed_gatk_files/input_1_pe/variantCounts_filtered_by_library_err_corrected_false_doubles.csv",
                                      aa_seq_file_path = "intermediate_files/aa_seq.txt",
                                      output_csv_path = "intermediate_files/processed_gatk_files/input_1_pe/variantCounts_for_heatmaps_err_corrected_false_doubles.csv")

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
