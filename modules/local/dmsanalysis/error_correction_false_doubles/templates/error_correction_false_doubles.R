#!/usr/bin/env Rscript

## false-doubles based sequencing error correction in nf-core/deepmutscan
## The library only contains single-codon changes, so observed multi-codon variants are
## sequencing errors; their counts are used to estimate and subtract the per-1nt error rate.

## --- Helper functions ---

# false-doubles based sequencing error correction (nucleotide level)
seq_error_correct_by_false_doubles <- function(input_count_path_raw, input_count_path_processed, output_file_path){

  # raw (pre-library-filter) counts have no header
  colnames <- c("counts", "cov", "mean_length_variant_reads", "varying_bases",
                "base_mut", "varying_codons", "codon_mut", "aa_mut", "pos_mut")
  input.counts.raw <- read.table(input_count_path_raw, sep = "\\t", header = FALSE, fill = TRUE, col.names = colnames)
  input.counts.processed <- read.csv(input_count_path_processed)

  ## append key columns
  input.counts.processed\$counts_corrected <- input.counts.processed\$counts
  input.counts.processed\$counts_per_cov_corrected <- input.counts.processed\$counts_per_cov

  cat("Sequencing error correction (false doubles)...\\n")
  for (i in grep("[,]", input.counts.processed[,"base_mut"], invert = T)){

    tmp.single <- input.counts.processed[i,"base_mut"]

    ## locate this variant across all multi-codon variants in the original count matrix
    tmp.false.doubles <- input.counts.raw[grep(paste0("(?<!\\\\d)",input.counts.processed[i,"codon_mut"]), input.counts.raw[,"codon_mut"], perl = T),]

    ## subset false multi-codon variants
    tmp.false.doubles <- tmp.false.doubles[which(tmp.false.doubles[,"varying_bases"] >= 3 & tmp.false.doubles[,"varying_codons"] < 3),]
    if(nrow(tmp.false.doubles) == 0){
      next
    }

    ## match with the corresponding correct single codon variant(s)
    tmp.true.singles <- tmp.false.doubles\$base_mut
    tmp.true.singles <- strsplit(tmp.true.singles, ", ")
    tmp.true.singles <- lapply(tmp.true.singles, function(x){x <- x[x != tmp.single]; x <- paste0(x, collapse = ", "); return(x)})
    tmp.true.singles <- do.call(c, tmp.true.singles)
    tmp.true.singles <- input.counts.processed[match(tmp.true.singles, input.counts.processed[,"base_mut"]),]
    if(all(is.na(tmp.true.singles) == T)){
      next
    }

    ## only count events for which there are both true single codon and false double codon variants
    if(any(is.na(tmp.true.singles\$counts))){
      tmp.false.doubles <- tmp.false.doubles[-which(is.na(tmp.true.singles\$counts)),]
      tmp.true.singles <- tmp.true.singles[-which(is.na(tmp.true.singles\$counts)),]
    }
    if(all(is.na(tmp.true.singles) == T)){
      next
    }

    ## expected false 1nt count probability (MLE)
    e_MLE <- sum(tmp.false.doubles\$counts) / sum(tmp.true.singles\$counts * c(tmp.false.doubles\$cov / tmp.true.singles\$cov))
    input.counts.processed[i,"counts_per_cov_corrected"] <- input.counts.processed[i,"counts_per_cov"] - e_MLE

    ## adjust the total counts by cross-multiplication
    input.counts.processed[i,"counts_corrected"] <- input.counts.processed[i,"counts"] * c(input.counts.processed[i,"counts_per_cov_corrected"] / input.counts.processed[i,"counts_per_cov"])

    ## round to nearest integer, do not allow for negative counts
    input.counts.processed[i,"counts_corrected"] <- round(input.counts.processed[i,"counts_corrected"])
    if(input.counts.processed[i,"counts_corrected"] < 0){
      input.counts.processed[i,"counts_corrected"] <- 0
      input.counts.processed[i,"counts_per_cov_corrected"] <- 0
    }

  }

  ## corrected counts become the canonical columns so every count-dependent step uses them;
  ## the originals are preserved as *_raw for the error-correction report
  input.counts.processed\$counts_raw <- input.counts.processed\$counts
  input.counts.processed\$counts_per_cov_raw <- input.counts.processed\$counts_per_cov
  input.counts.processed\$counts <- input.counts.processed\$counts_corrected
  input.counts.processed\$counts_per_cov <- input.counts.processed\$counts_per_cov_corrected

  write.csv(input.counts.processed, file = output_file_path, row.names = F)

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
seq_error_correct_by_false_doubles(
  input_count_path_raw = "$raw_counts",
  input_count_path_processed = "$filtered",
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
