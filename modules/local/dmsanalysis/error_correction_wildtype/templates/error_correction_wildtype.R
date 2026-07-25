#!/usr/bin/env Rscript

## Wildtype-based sequencing error correction in nf-core/deepmutscan
## An additional deep wildtype-only sequencing sample measures the position-specific error
## profile, which is subtracted from each library sample's per-base coverage.

## --- Helper functions ---

# WT-based sequencing error correction (nucleotide level)
seq_error_correct_by_WT_nt <- function(wt_seq_count_path, input_count_path, output_file_path){

  WT.counts <- read.csv(wt_seq_count_path)
  input.counts <- read.csv(input_count_path)

  ## append key columns
  input.counts\$counts_corrected <- input.counts\$counts
  input.counts\$counts_per_cov_corrected <- input.counts\$counts_per_cov

  for (i in 1:nrow(WT.counts)){

    ## only look at variants observed in both the input and WT sequencing
    tmp.id <- match(WT.counts[i,"base_mut"], input.counts[,"base_mut"])
    if(is.na(tmp.id) == T){
      next
    }

    ## subtract the observed per-base coverage in the WT sequencing
    input.counts[tmp.id,"counts_per_cov_corrected"] <- input.counts[tmp.id,"counts_per_cov"] - WT.counts[i,"counts_per_cov"]

    ## adjust the total counts by cross-multiplication
    input.counts[tmp.id,"counts_corrected"] <- input.counts[tmp.id,"counts"] * c(input.counts[tmp.id,"counts_per_cov_corrected"] / input.counts[tmp.id,"counts_per_cov"])

    ## round to nearest integer, do not allow for negative counts
    input.counts[tmp.id,"counts_corrected"] <- round(input.counts[tmp.id,"counts_corrected"])
    if(input.counts[tmp.id,"counts_corrected"] < 0){
      input.counts[tmp.id,"counts_corrected"] <- 0
      input.counts[tmp.id,"counts_per_cov_corrected"] <- 0
    }

  }

  ## corrected counts become the canonical columns so every count-dependent step uses them;
  ## the originals are preserved as *_raw for the error-correction report
  input.counts\$counts_raw <- input.counts\$counts
  input.counts\$counts_per_cov_raw <- input.counts\$counts_per_cov
  input.counts\$counts <- input.counts\$counts_corrected
  input.counts\$counts_per_cov <- input.counts\$counts_per_cov_corrected

  write.csv(input.counts, file = output_file_path, row.names = F)

}

# propagate corrected counts onto the completed library table (used by e.g. logdiff)
correct_library_completed <- function(completed_path, corrected_filtered_path, output_file_path){
  completed <- read.csv(completed_path)
  corrected <- read.csv(corrected_filtered_path)
  idx <- match(completed\$codon_mut, corrected\$codon_mut)
  has <- !is.na(idx)
  completed\$counts[has] <- corrected\$counts[idx[has]]
  completed\$counts_per_cov[has] <- corrected\$counts_per_cov[idx[has]]
  write.csv(completed, file = output_file_path, row.names = FALSE)
}

# re-make the AA-level heatmap table from corrected counts (canonical column names)
seq_error_correct_counts_for_heatmaps <- function(input_csv_path, aa_seq_file_path, output_csv_path, threshold = 3) {

  raw_counts <- read.table(input_csv_path, sep = ",", header = TRUE, stringsAsFactors = FALSE)

  wt_seq <- readLines(aa_seq_file_path)
  wt_seq <- unlist(strsplit(wt_seq, ""))
  wt_seq[wt_seq == "X"] <- "*"

  # summarise corrected counts for each unique pos_mut
  aggregated_counts_per_cov <- aggregate(
    counts_per_cov_corrected ~ pos_mut,
    data = raw_counts,
    FUN = function(x) sum(x, na.rm = TRUE)
  )
  aggregated_counts <- aggregate(
    counts_corrected ~ pos_mut,
    data = raw_counts,
    FUN = function(x) sum(x, na.rm = TRUE)
  )
  aggregated_data <- merge(aggregated_counts_per_cov, aggregated_counts, by = "pos_mut", all = TRUE)

  # canonical output names so the count heatmaps consume this drop-in
  names(aggregated_data)[names(aggregated_data) == "counts_per_cov_corrected"] <- "total_counts_per_cov"
  names(aggregated_data)[names(aggregated_data) == "counts_corrected"] <- "total_counts"

  aggregated_data\$wt_aa <- sub("(\\\\D)(\\\\d+)(\\\\D)", "\\\\1", aggregated_data\$pos_mut)
  aggregated_data\$position <- as.numeric(sub("(\\\\D)(\\\\d+)(\\\\D)", "\\\\2", aggregated_data\$pos_mut))
  aggregated_data\$mut_aa <- sub("(\\\\D)(\\\\d+)(\\\\D)", "\\\\3", aggregated_data\$pos_mut)
  aggregated_data\$mut_aa[aggregated_data\$mut_aa == "X"] <- "*"

  all_amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                       "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "*")
  all_positions <- seq_along(wt_seq)
  complete_data <- expand.grid(mut_aa = all_amino_acids, position = all_positions, stringsAsFactors = FALSE)

  heatmap_data <- merge(complete_data, aggregated_data, by = c("mut_aa", "position"), all.x = TRUE, sort = FALSE)

  heatmap_data\$total_counts_per_cov[is.na(heatmap_data\$total_counts_per_cov)] <- 0
  heatmap_data\$wt_aa <- wt_seq[heatmap_data\$position]

  low_count <- is.na(heatmap_data\$total_counts) | heatmap_data\$total_counts < threshold
  heatmap_data\$total_counts_per_cov[low_count] <- NA
  heatmap_data\$total_counts[low_count] <- NA

  missing_pos_mut <- is.na(heatmap_data\$pos_mut)
  heatmap_data\$pos_mut[missing_pos_mut] <- paste0(
    heatmap_data\$wt_aa[missing_pos_mut],
    heatmap_data\$position[missing_pos_mut],
    heatmap_data\$mut_aa[missing_pos_mut]
  )

  out.order <- paste0(rep(1:max(heatmap_data\$position), each = 21),
                      rep(c("A", "C", "D", "E", "F", "G", "H",
                            "I", "K", "L","M", "N", "P", "Q",
                            "R", "S", "T", "V", "W", "Y", "*"),
                          max(heatmap_data\$position)))
  heatmap_data <- heatmap_data[match(out.order, paste0(heatmap_data\$position, heatmap_data\$mut_aa)),]
  rownames(heatmap_data) <- 1:nrow(heatmap_data)

  write.csv(heatmap_data, file = output_csv_path, row.names = FALSE)

}

## --- Run ---
seq_error_correct_by_WT_nt(
  wt_seq_count_path = "$wt_counts",
  input_count_path = "$filtered",
  output_file_path = "variantCounts_filtered_by_library_error_corrected.csv"
)

seq_error_correct_counts_for_heatmaps(
  input_csv_path = "variantCounts_filtered_by_library_error_corrected.csv",
  aa_seq_file_path = "$aa_seq",
  output_csv_path = "variantCounts_for_heatmaps_error_corrected.csv",
  threshold = $min_counts
)

correct_library_completed(
  completed_path = "$completed",
  corrected_filtered_path = "variantCounts_filtered_by_library_error_corrected.csv",
  output_file_path = "library_completed_variantCounts_error_corrected.csv"
)

## --- versions ---
r_version <- strsplit(version[['version.string']], ' ')[[1]][3]
if (is.null(r_version)) r_version <- "unknown"
f <- file("versions.yml", "w")
writeLines(
  c(
    '"${task.process}":',
    paste('    r-base:', r_version)
  ),
  f
)
close(f)
