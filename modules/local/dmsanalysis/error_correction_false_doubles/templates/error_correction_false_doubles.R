#!/usr/bin/env Rscript
## nf-core/deepmutscan - false-double-mutant (FDM) sequencing-error correction, one library.
##
## Applies the per-SNV error rates of the library's pool: e_j * coverage_j of this library is subtracted
## from the SNV's count (rounded, floored at zero). SNVs without a rate are left as they are;
## multi-nucleotide codon variants are not corrected. The reference must be in the frame of the counts
## (checked first). The correction is the fdmcorrect R package (bioconda::dms-fdmcorrect,
## https://github.com/MaximilianStammnitz/dms-fdmcorrect); the corrected counts become the canonical
## columns, with the originals kept as *_raw for the error-correction report.
##
## deepmutscan's own two steps follow: the amino-acid table for the count heatmaps, and the corrected
## counts on the completed library table.
##
## Nextflow template: an R `\$` is a literal dollar; backslashes meant for R are doubled.

suppressPackageStartupMessages(library(fdmcorrect))

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

  # summarise corrected counts for each unique pos_mut (a library without any programmed variant gives
  # an empty table, and every cell of the heatmap table is then missing)
  if (nrow(raw_counts) > 0) {
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
  } else {
    aggregated_data <- data.frame(pos_mut = character(0), counts_per_cov_corrected = numeric(0),
                                  counts_corrected = numeric(0), stringsAsFactors = FALSE)
  }

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

## -------------------------------------------------------------------------------------------------
## Run
## -------------------------------------------------------------------------------------------------
fdm_correct(
    library_variants = "$filtered",
    evidence         = "$evidence",
    rates            = "$rates",
    reference        = "$wt_fasta",
    outdir           = ".",
    out              = "variantCounts_filtered_by_library_error_corrected.csv",
    seq_error_rate   = "seq_error_rate.csv"
)

seq_error_correct_counts_for_heatmaps(
  input_csv_path = "variantCounts_filtered_by_library_error_corrected.csv",
  aa_seq_file_path = "$aa_seq",
  output_csv_path = "variantCounts_for_heatmaps_error_corrected.csv",
  threshold = 1  # every observed variant; fitness applies --min_counts to the input libraries
)

correct_library_completed(
  completed_path = "$completed",
  corrected_filtered_path = "variantCounts_filtered_by_library_error_corrected.csv",
  output_file_path = "library_completed_variantCounts_error_corrected.csv"
)

## --- versions ---
r_version <- strsplit(version[['version.string']], ' ')[[1]][3]
if (is.null(r_version)) r_version <- "unknown"
writeLines(c('"${task.process}":', paste("    r-base:", r_version),
             paste("    dms-fdmcorrect:", as.character(packageVersion("fdmcorrect")))), "versions.yml")
