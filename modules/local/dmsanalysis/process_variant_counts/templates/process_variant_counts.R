#!/usr/bin/env Rscript

# four parts condensed into one R script

#########
# process raw counts
#########

# input: counts variantcounts tsv file, output_path
# output: csv with column names. Creates additional counts_per_cov column. Fills pos_mut column for synonymous mutations. Sorted out variants that have mutations, but do not show up in the specifying columns -> was affecting roughly 30 low-count variants out of over 15000 in Taylor's data

process_raw_counts <- function(counts_file_path, output_csv_path) {

  # Set the column names
  colnames <- c("counts", "cov", "mean_length_variant_reads", "varying_bases",
                          "base_mut", "varying_codons", "codon_mut", "aa_mut", "pos_mut")

  # Read the variant-count file into a data frame
  counts_raw <- read.table(counts_file_path, sep = "\\t", header = FALSE, fill = TRUE, col.names = colnames)

  # Filter out rows where 'aa_mut' is empty or NA
  counts_raw <- counts_raw[!(counts_raw\$aa_mut == "" | is.na(counts_raw\$aa_mut)), ]

  # Handle synonymous mutations: where aa_mut starts with "S" and pos_mut is either NA or ""
  # (ifelse, paste0 and sub work on whole columns, so no row-by-row evaluation is needed)
  counts_raw\$pos_mut <- ifelse(
    (is.na(counts_raw\$pos_mut) | counts_raw\$pos_mut == "") & grepl("^S:", counts_raw\$aa_mut),
    # Construct the new 'pos_mut' entry for synonymous mutations
    paste0(
      sub("S:([A-Z])>[A-Z]", "\\\\1", counts_raw\$aa_mut),  # Get the original amino acid from 'aa_mut'
      sub("^(\\\\d+):.*", "\\\\1", counts_raw\$codon_mut),    # Get the position from 'codon_mut'
      sub("S:[A-Z]>([A-Z])", "\\\\1", counts_raw\$aa_mut)   # Get the mutated amino acid from 'aa_mut'
    ),
    counts_raw\$pos_mut  # Keep the existing 'pos_mut' if it's not NA or ""
  )
  counts_raw\$counts_per_cov <- counts_raw\$counts / counts_raw\$cov

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

filter_counts_by_codon_library <- function(counts_file_path, codon_library_path, output_file_path) {
  # Load the variant-count table from the provided file path
  counts_table <- read.csv(counts_file_path)

  # Load the codon library from the provided .csv file
  codon_library <- read.csv(codon_library_path)

  # Ensure the codon library has the expected columns
  if (!all(c("Codon_Number", "wt_codon", "Variant") %in% colnames(codon_library))) {
    stop("Codon library must contain columns 'Codon_Number', 'wt_codon', and 'Variant'.")
  }

  # Filter the variant-count table: single-codon variants whose position and mutated codon are in
  # the library (as a programmed variant or as the wildtype codon), and whose changed bases lie
  # within one codon (first to last changed position at most 2 apart). Both checks are evaluated
  # for whole columns at once; row-by-row evaluation (rowwise()) was the slow part for long ORFs.
  counts_table\$codon_mut <- toupper(counts_table\$codon_mut)
  library_keys <- unique(c(
    paste(codon_library\$Codon_Number, codon_library\$Variant),
    paste(codon_library\$Codon_Number, codon_library\$wt_codon)
  ))
  # the first match of a regular expression in each string, NA where there is none (stringr::str_extract)
  extract_first_match <- function(x, pattern) {
    m <- regexpr(pattern, x)
    out <- rep(NA_character_, length(x))
    hit <- !is.na(m) & m != -1
    out[hit] <- substr(x[hit], m[hit], m[hit] + attr(m, "match.length")[hit] - 1L)
    out
  }
  changed_base_span <- function(base_mut) {
    vapply(strsplit(base_mut, ",\\\\s*"), function(mutations) {
      positions <- suppressWarnings(as.numeric(extract_first_match(mutations, "^[0-9]+")))
      suppressWarnings(max(positions, na.rm = TRUE) - min(positions, na.rm = TRUE))
    }, numeric(1))
  }
  # (each step keeps the rows where its condition is TRUE, as dplyr::filter() did: NA drops the row)
  filtered_counts <- counts_table[which(counts_table\$varying_codons == 1), , drop = FALSE]
  filtered_counts <- filtered_counts[which(paste(as.numeric(sub(":.*", "", filtered_counts\$codon_mut)),
                                                 sub(".*>", "", filtered_counts\$codon_mut)) %in% library_keys), , drop = FALSE]
  filtered_counts <- filtered_counts[which(changed_base_span(filtered_counts\$base_mut) <= 2), , drop = FALSE]

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
# output: completed counts_file with all possible variants, including those never measured in sequencing (counts and counts_per_cov 0, observed FALSE)

# Row indices of a join on one key column (dplyr 1.1 semantics): every row of x in order, each with
# all matching rows of y in y order (NA where there is none); NA keys match each other.
join_rows <- function(x_key, y_key) {
  u <- unique(y_key)
  y_by_key <- split(seq_along(y_key), match(y_key, u))
  g <- match(x_key, u)
  n <- ifelse(is.na(g), 1L, lengths(y_by_key)[g])
  x_rows <- rep(seq_along(x_key), n)
  y_rows <- unlist(lapply(seq_along(g), function(i) if (is.na(g[i])) NA_integer_ else y_by_key[[g[i]]]),
                   use.names = FALSE)
  if (is.null(y_rows)) y_rows <- integer(0)
  list(x = x_rows, y = y_rows)
}

# The columns of x (rows x_rows) followed by the columns of y except `by` (rows y_rows); NA rows give
# NA of the column's type
bind_join <- function(x, y, x_rows, y_rows, by) {
  out <- lapply(x, function(col) col[x_rows])
  for (col in setdiff(names(y), by)) out[[col]] <- y[[col]][y_rows]
  as.data.frame(out, stringsAsFactors = FALSE, optional = TRUE)
}

# Hamming distance between two codons (varying_bases), for whole columns at once
hamming_distance <- function(wt_codon, variant_codon) {
  Reduce(`+`, lapply(1:3, function(i) as.integer(substr(wt_codon, i, i) != substr(variant_codon, i, i))))
}

# The standard genetic code (Biostrings::GENETIC_CODE); stop codons are "*"
GENETIC_CODE_STANDARD <- c(
  TTT = "F", TTC = "F", TTA = "L", TTG = "L", CTT = "L", CTC = "L", CTA = "L", CTG = "L",
  ATT = "I", ATC = "I", ATA = "I", ATG = "M", GTT = "V", GTC = "V", GTA = "V", GTG = "V",
  TCT = "S", TCC = "S", TCA = "S", TCG = "S", CCT = "P", CCC = "P", CCA = "P", CCG = "P",
  ACT = "T", ACC = "T", ACA = "T", ACG = "T", GCT = "A", GCC = "A", GCA = "A", GCG = "A",
  TAT = "Y", TAC = "Y", TAA = "*", TAG = "*", CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
  AAT = "N", AAC = "N", AAA = "K", AAG = "K", GAT = "D", GAC = "D", GAA = "E", GAG = "E",
  TGT = "C", TGC = "C", TGA = "*", TGG = "W", CGT = "R", CGC = "R", CGA = "R", CGG = "R",
  AGT = "S", AGC = "S", AGA = "R", AGG = "R", GGT = "G", GGC = "G", GGA = "G", GGG = "G"
)

# Amino acid of each codon (NA for invalid codons), for whole columns at once
get_amino_acid <- function(codon) {
  unname(GENETIC_CODE_STANDARD[toupper(codon)])
}

complete_prefiltered_counts <- function(possible_nnk_path, prefiltered_counts_path, output_file_path) {

  # Load the possible NNK mutations CSV
  possible_nnk <- read.csv(possible_nnk_path)

  # Load the prefiltered variant-count CSV
  prefiltered_counts <- read.csv(prefiltered_counts_path)

  # Create codon_mut column in possible_NNK_mutations in the format 'Codon_Number:wt_codon>Variant'
  possible_nnk\$codon_mut <- paste0(possible_nnk\$Codon_Number, ":", possible_nnk\$wt_codon, ">", possible_nnk\$Variant)

  # Merge both dataframes based on the codon_mut column (full join to include all): every counted
  # variant with its library row, then the library rows that were never counted
  prefiltered_counts\$codon_mut <- toupper(prefiltered_counts\$codon_mut)
  j <- join_rows(prefiltered_counts\$codon_mut, possible_nnk\$codon_mut)
  library_only <- setdiff(seq_len(nrow(possible_nnk)), j\$y)
  merged_data <- bind_join(prefiltered_counts, possible_nnk,
                           c(j\$x, rep(NA_integer_, length(library_only))), c(j\$y, library_only), "codon_mut")
  merged_data\$codon_mut[length(j\$x) + seq_along(library_only)] <- possible_nnk\$codon_mut[library_only]

  # Programmed variants that were never observed get their true count of 0, flagged with
  # observed = FALSE. Plots on a log scale (the LogDiff plot) place zeros themselves, at a floor
  # derived from the data; a fixed placeholder (previously 1e-6) can rank above real variants in
  # deeply sequenced libraries, where one read is fewer than 1e-6 counts per coverage.
  merged_data\$observed <- !is.na(merged_data\$counts)
  merged_data\$counts_per_cov <- ifelse(is.na(merged_data\$counts_per_cov), 0, merged_data\$counts_per_cov)
  merged_data\$counts <- ifelse(is.na(merged_data\$counts), 0, merged_data\$counts)

  # Hamming distance (varying_bases) and mutation details: aa_mut (M = missense, S = synonymous,
  # e.g. M:D>S) and pos_mut (wildtype AA, codon position, mutated AA, e.g. D2Q)
  wt_aa      <- get_amino_acid(merged_data\$wt_codon)
  variant_aa <- get_amino_acid(merged_data\$Variant)
  merged_data\$varying_bases <- hamming_distance(merged_data\$wt_codon, merged_data\$Variant)
  merged_data\$aa_mut <- paste0(ifelse(wt_aa != variant_aa, "M", "S"), ":", wt_aa, ">", variant_aa)
  merged_data\$pos_mut <- paste0(wt_aa, merged_data\$Codon_Number, variant_aa)

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

prepare_counts_data_for_counts_heatmaps <- function(counts_file_path, aa_seq_file_path, output_csv_path, threshold = 3) {
  # Load the raw variant-count data
  raw_counts <- read.table(counts_file_path, sep = ",", header = TRUE)

  # Read the wild-type amino acid sequence from the text file
  wt_seq <- readLines(aa_seq_file_path)
  wt_seq <- unlist(strsplit(wt_seq, ""))  # Split the sequence into individual amino acids

  # Summarize counts-per-cov for each unique aa mutation in pos_mut: one row per pos_mut, sorted in
  # the C locale with NA last (as dplyr's group_by() sorts), sums over the rows in table order
  groups <- unique(raw_counts\$pos_mut)
  groups <- groups[order(groups, method = "radix", na.last = TRUE)]
  rows_by_group <- unname(split(seq_len(nrow(raw_counts)), factor(match(raw_counts\$pos_mut, groups), levels = seq_along(groups))))
  aggregated_data <- data.frame(
    pos_mut = groups,
    total_counts_per_cov = vapply(rows_by_group, function(i) sum(raw_counts\$counts_per_cov[i], na.rm = TRUE),
                                  sum(raw_counts\$counts_per_cov[0])),
    total_counts = vapply(rows_by_group, function(i) sum(raw_counts\$counts[i], na.rm = TRUE),  # Also sum the counts
                          sum(raw_counts\$counts[0])),
    stringsAsFactors = FALSE
  )

  # Extract the wild-type position and mutations from 'pos_mut'
  aggregated_data\$wt_aa <- sub("(\\\\D)(\\\\d+)(\\\\D)", "\\\\1", aggregated_data\$pos_mut)  # Wild-type amino acid (e.g., S)
  aggregated_data\$position <- as.numeric(sub("(\\\\D)(\\\\d+)(\\\\D)", "\\\\2", aggregated_data\$pos_mut))  # Position (e.g., 3)
  aggregated_data\$mut_aa <- sub("(\\\\D)(\\\\d+)(\\\\D)", "\\\\3", aggregated_data\$pos_mut)   # Mutant amino acid (e.g., R)

  # Replace 'X' with '*', indicating the stop codon
  aggregated_data\$mut_aa <- ifelse(aggregated_data\$mut_aa == "X", "*", aggregated_data\$mut_aa)

  # Replace 'X' with '*' in the wild-type amino acid sequence as well
  wt_seq <- ifelse(wt_seq == "X", "*", wt_seq)

  # Define all 20 standard amino acids and the stop codon "*"
  all_amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                       "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "*")

  # Create a list of all positions in the wild-type sequence
  all_positions <- 1:length(wt_seq)

  # Create a complete grid of all possible combinations of positions and amino acids
  complete_data <- expand.grid(mut_aa = all_amino_acids, position = all_positions,
                               KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  complete_data\$position <- as.numeric(complete_data\$position)  # the joined column is numeric

  # Merge the summarized data with the complete grid (filling missing entries with 0): every grid cell
  # in grid order, with its row of the summary (by mutant amino acid and position)
  key_sep <- "\\r"
  j <- join_rows(paste(complete_data\$mut_aa, complete_data\$position, sep = key_sep),
                 paste(aggregated_data\$mut_aa, aggregated_data\$position, sep = key_sep))
  heatmap_data <- bind_join(complete_data, aggregated_data[, setdiff(names(aggregated_data), c("mut_aa", "position")), drop = FALSE],
                            j\$x, j\$y, character(0))
  heatmap_data\$total_counts_per_cov <- ifelse(is.na(heatmap_data\$total_counts_per_cov), 0, heatmap_data\$total_counts_per_cov)
  heatmap_data\$wt_aa <- wt_seq[heatmap_data\$position]  # Assign the wild-type amino acid

  # Set variants with counts < threshold to NA
  heatmap_data\$total_counts_per_cov <- ifelse(heatmap_data\$total_counts < threshold, NA, heatmap_data\$total_counts_per_cov)
  heatmap_data\$total_counts <- ifelse(heatmap_data\$total_counts < threshold, NA, heatmap_data\$total_counts)

  # Fill pos_mut column
  heatmap_data\$pos_mut <- ifelse(is.na(heatmap_data\$pos_mut),
                                 paste0(heatmap_data\$wt_aa, heatmap_data\$position, heatmap_data\$mut_aa),
                                 heatmap_data\$pos_mut)

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
  threshold = 1  # every observed variant; fitness applies --min_counts to the input libraries
)







####
# create versions.yml
####
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
