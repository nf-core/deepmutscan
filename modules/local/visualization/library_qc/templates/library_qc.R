#!/usr/bin/env Rscript

# Library QC of one library in one task: the six plot sets that used to be six processes
# (VISUALIZATION_COUNTS_HEATMAP, _COUNTS_PER_COV, _GLOBAL_POS_BIASES_COUNTS, _GLOBAL_POS_BIASES_COV,
# _LOGDIFF and _SEQDEPTH), run in sequence in one R session. Each plot function below is that
# process's script; the two sliding-window functions were renamed (both were called `position_biases`),
# and the dplyr and zoo calls are base R (helpers below, with the same results). Every output file is
# the same as before.
#
# Inputs:  variantCounts_filtered_by_library (per-variant counts, coverage and counts per coverage
#          of the programmed variants), variantCounts_for_heatmaps (amino-acid level table),
#          library_completed_variantCounts (every programmed variant, with 0 for the unobserved),
#          aa_seq (wildtype amino-acid sequence), possible_mutations (the programmed codon variants,
#          only read for the rarefaction), sliding_window_size, aimed_cov, run_seqdepth
# Outputs: counts_heatmap.pdf, counts_per_cov_heatmap.pdf, rolling_counts.pdf,
#          rolling_counts_per_cov.pdf, rolling_coverage.pdf, logdiff_plot.pdf,
#          logdiff_varying_bases.pdf and - with run_seqdepth - SeqDepth.pdf and seqdepth_curve.csv

library(ggplot2)
library(scales)

####################################################################################################
# Base-R helpers (for what dplyr and zoo did before, with the same results)
####################################################################################################

# dplyr::bind_rows(a, b): the columns of both in order of appearance, missing ones as NA; a factor
# combined with a character column becomes character
bind_rows_base <- function(a, b) {
  cols <- union(names(a), names(b))
  out <- lapply(cols, function(cn) {
    va <- a[[cn]]; vb <- b[[cn]]
    if (is.factor(va) && is.character(vb)) va <- as.character(va)
    if (is.factor(vb) && is.character(va)) vb <- as.character(vb)
    if (is.null(va)) va <- vb[rep(NA_integer_, nrow(a))]
    if (is.null(vb)) vb <- va[rep(NA_integer_, nrow(b))]
    c(va, vb)
  })
  names(out) <- cols
  as.data.frame(out, stringsAsFactors = FALSE, optional = TRUE)
}

# The rows of d in a stable ascending order of x, NA last (dplyr::arrange)
arrange_base <- function(d, x) {
  d <- d[order(x, method = "radix"), , drop = FALSE]
  rownames(d) <- NULL
  d
}

# The rows of d where cond is TRUE, NA dropping the row (dplyr::filter)
filter_base <- function(d, cond) {
  d <- d[which(cond), , drop = FALSE]
  rownames(d) <- NULL
  d
}

# zoo::rollapply(x, width, FUN, fill = "extend"): FUN over each centred window that lies fully
# inside x (offsets floor(width/2) - width + 1 ... floor(width/2)), then zoo::na.fill(, "extend"):
# the first result repeated before, the last after, and inside the range a linear interpolation
# over the windows whose result is not NA (stats::approx), exactly as zoo 1.8 does it
roll_extend <- function(x, width, FUN) {
  n <- length(x)
  if (n < 1L) return(x)
  w <- trunc(width)
  offsets <- seq(to = floor(w / 2), length.out = w)
  res <- vector("list", n)
  for (i in seq_len(n)) {
    p <- i + offsets
    if (all(p >= 1 & p <= n)) res[[i]] <- FUN(x[p])
  }
  ix <- !vapply(res, is.null, logical(1))
  out <- do.call("c", lapply(res, function(r) if (is.null(r)) NA else r))
  if (any(ix)) {
    wix <- which(ix)
    lo <- wix[1]
    hi <- wix[length(wix)]
    rng <- seq(lo, length.out = hi - lo + 1)
    out[seq_len(lo - 1)] <- out[[lo]]
    out[rng] <- approx(wix, out[wix], xout = rng)\$y
    out[seq(hi + 1, length.out = n - hi)] <- out[[hi]]
  } else {
    out[is.na(out)] <- NA
  }
  out
}

####################################################################################################
# 1. Heatmap of counts per variant (counts_heatmap.pdf)
####################################################################################################
# Input: variantCounts_for_heatmaps table, output path, dropout threshold
# Output: counts_heatmap.pdf

counts_heatmap <- function(input_csv_path, threshold = 3, output_pdf_path, img_format = "pdf") {

  # Inner function to add padding to the last row, adding 21 amino acids per position
  pad_heatmap_data_long <- function(heatmap_data_long, min_non_na_value, num_positions_per_row = 75) {
    all_amino_acids <- c("G", "A", "V", "L", "M", "I", "F",
                         "Y", "W", "K", "R", "H", "D", "E",
                         "S", "T", "C", "N", "Q", "P", "*")

    max_position <- max(heatmap_data_long\$position)
    num_missing_positions <- num_positions_per_row - (max_position %% num_positions_per_row)

    if (num_missing_positions < num_positions_per_row) {
      new_positions <- (max_position + 1):(max_position + num_missing_positions)

      # Add all 21 amino acid variants for each new position
      padding_data <- expand.grid(
        mut_aa = all_amino_acids,  # All possible amino acids
        position = new_positions  # New positions to be padded
      )

      # Set placeholder values for the added positions to the exact smallest non-NA value
      padding_data\$total_counts <- min_non_na_value  # Set to the smallest non-NA value
      padding_data\$wt_aa <- "Y"  # Set wild-type amino acid to 'Y'
      padding_data\$wt_aa_pos <- paste0("Y", padding_data\$position)  # Create wt_aa_pos with correct positions
      padding_data\$row_group <- max(heatmap_data_long\$row_group)  # Set row group to the current last group

      # Add the new padding rows to heatmap_data_long
      heatmap_data_long <- bind_rows_base(heatmap_data_long, padding_data)
    }

    return(heatmap_data_long)
  }

  # Load the CSV data
  heatmap_data <- read.csv(input_csv_path)

  # Check if the necessary column exists in the data
  if (!"total_counts" %in% colnames(heatmap_data)) {
    stop("The column 'total_counts' is not found in the data.")
  }

  # Create heatmap_data_long by selecting necessary columns
  heatmap_data_long <- heatmap_data[, c("mut_aa", "position", "total_counts", "wt_aa")]  # Use 'total_counts'

  # Find the smallest non-NA value in total_counts
  min_non_na_value <- min(heatmap_data_long\$total_counts, na.rm = TRUE)

  # Group positions by rows (75 positions per row) and calculate row_group
  heatmap_data_long\$row_group <- ((heatmap_data_long\$position - 1) %/% 75) + 1  # Grouping positions into rows

  # Apply padding to add missing positions at the end of the last row, using the calculated min value
  heatmap_data_long <- pad_heatmap_data_long(heatmap_data_long, min_non_na_value)

  # Convert positions to numeric, sort them, and create wt_aa_pos for the plot
  heatmap_data_long\$position <- as.numeric(heatmap_data_long\$position)  # Ensure position is numeric
  heatmap_data_long <- arrange_base(heatmap_data_long, heatmap_data_long\$position)  # Sort by position
  heatmap_data_long\$wt_aa_pos <- factor(paste0(heatmap_data_long\$wt_aa, heatmap_data_long\$position),
                                       levels = unique(paste0(heatmap_data_long\$wt_aa, heatmap_data_long\$position)))  # Create sorted factor levels for wt_aa_pos

  # Add a column to identify synonymous mutations (where mut_aa == wt_aa)
  heatmap_data_long\$synonymous <- heatmap_data_long\$mut_aa == heatmap_data_long\$wt_aa

  # Definiere die korrekte Reihenfolge der Aminosäuren
  amino_acid_order <- rev(c("G", "A", "V", "L", "M", "I", "F",
                         "Y", "W", "K", "R", "H", "D", "E",
                         "S", "T", "C", "N", "Q", "P", "*"))

  heatmap_data_long\$mut_aa <- factor(heatmap_data_long\$mut_aa, levels = amino_acid_order)

  # Bearbeite heatmap_data_long und erstelle syn_positions gleichzeitig
  syn_positions <- heatmap_data_long
  syn_positions\$mut_aa <- factor(syn_positions\$mut_aa, levels = amino_acid_order)
  # Berechne die x-Koordinate, die pro Gruppe immer von 1 bis 75 verläuft
  syn_positions\$x <- as.numeric(factor(syn_positions\$wt_aa_pos, levels = unique(syn_positions\$wt_aa_pos))) - ((syn_positions\$row_group - 1) * 75)
  syn_positions\$y <- as.numeric(factor(syn_positions\$mut_aa, levels = amino_acid_order))
  syn_positions <- filter_base(syn_positions, syn_positions\$synonymous == TRUE)

  # Calculate the number of row groups and adjust plot height dynamically
  num_row_groups <- max(heatmap_data_long\$row_group)
  plot_height <- num_row_groups * 4

  # Set the limits for the color scale, ignoring NA (negative values are now NA)
  min_count <- min(heatmap_data_long\$total_counts, na.rm = TRUE)
  max_count <- max(heatmap_data_long\$total_counts, na.rm = TRUE)
  max_position <- max(heatmap_data\$position)

  # Create the heatmap plot with explicit handling for positions > max_position
  heatmap_plot <- ggplot(heatmap_data_long, aes(x = wt_aa_pos, y = mut_aa, fill = total_counts)) +
    scale_fill_gradientn(colours = c(alpha("blue", 0), "blue"), na.value = "grey35", trans = "log",  # Apply log transformation to the scale
                         limits = c(min_count, max_count),
                         breaks = scales::trans_breaks("log10", function(x) 10^x),  # Logarithmic scale breaks
                         labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_x_discrete(labels = function(x) {
      numeric_pos <- as.numeric(gsub("[^0-9]", "", x))
      ifelse(numeric_pos > max_position, " ", x)
    }) +
    geom_tile() +

    # Add diagonal lines for synonymous mutations using geom_segment
    geom_segment(data = syn_positions[syn_positions\$position <= max_position, ],
                 aes(x = x - 0.485, xend = x + 0.485,
                     y = y - 0.485, yend = y + 0.485, color = synonymous),
                 size = 0.2) +

    # Manuelle Farbskala für die diagonalen Linien
    scale_color_manual(values = c("TRUE" = "grey10"), labels = c("TRUE" = "")) +

    theme_minimal() +
    labs(title = "Heatmap of Counts per Variant", x = "Wild-type Amino Acid", y = "Mutant Amino Acid", fill = "Counts", color = "Synonymous Mutation") +
    theme(plot.title = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),  # Larger y-axis labels
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.title = element_text(size = 14),  # Larger legend title
          legend.text = element_text(size = 12),  # Larger legend text
          panel.spacing = unit(0.1, "lines"),  # Adjust panel spacing
          strip.text = element_blank(),  # Remove row group labels (facet numbers)
          strip.background = element_blank(),
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank()) +  # Remove minor grid lines
    facet_wrap(~ row_group, scales = "free_x", ncol = 1) +  # Group by 75 positions per row
    theme(panel.spacing = unit(0.2, "lines"))

  heatmap_plot <- heatmap_plot +
    geom_point(data = heatmap_data_long, aes(size = ""), colour = "black", alpha = 0)  # Invisible points for legend
  heatmap_plot <- heatmap_plot +
    guides(size = guide_legend(paste("Dropout (Counts <", threshold, ")"), override.aes = list(shape = 15, size = 8, colour = "grey35", alpha = 1)))  # Define Legend for Dropouts

  # Save the heatmap plot
  if (img_format == "pdf") {
    ggsave(output_pdf_path, plot = heatmap_plot, width = 16, height = plot_height, dpi = 150, device = cairo_pdf)
  } else {
    ggsave(output_pdf_path, plot = heatmap_plot, width = 16, height = plot_height, dpi = 150)
  }

  if (file.exists(output_pdf_path)) {
    print("Heatmap image successfully created!")
  } else {
    print("Error: Heatmap image was not created.")
  }
}

####################################################################################################
# 2. Heatmap of counts per coverage (counts_per_cov_heatmap.pdf)
####################################################################################################
# Input: variantCounts_for_heatmaps table, output path, dropout threshold
# Output: counts_per_cov_heatmap.pdf

counts_per_cov_heatmap <- function(input_csv_path, threshold = 3, output_pdf_path, img_format = "pdf") {

  # Inner function to add padding to the last row, adding 21 amino acids per position
  pad_heatmap_data_long <- function(heatmap_data_long, min_non_na_value, num_positions_per_row = 75) {
    all_amino_acids <- c("G", "A", "V", "L", "M", "I", "F",
                         "Y", "W", "K", "R", "H", "D", "E",
                         "S", "T", "C", "N", "Q", "P", "*")

    max_position <- max(heatmap_data_long\$position)
    num_missing_positions <- num_positions_per_row - (max_position %% num_positions_per_row)

    if (num_missing_positions < num_positions_per_row) {
      new_positions <- (max_position + 1):(max_position + num_missing_positions)

      # Add all 21 amino acid variants for each new position
      padding_data <- expand.grid(
        mut_aa = all_amino_acids,  # All possible amino acids
        position = new_positions  # New positions to be padded
      )

      # Set placeholder values for the added positions to the exact smallest non-NA value
      padding_data\$total_counts_per_cov <- min_non_na_value  # Set to the smallest non-NA value
      padding_data\$wt_aa <- "Y"  # Set wild-type amino acid to 'Y'
      padding_data\$wt_aa_pos <- paste0("Y", padding_data\$position)  # Create wt_aa_pos with correct positions
      padding_data\$row_group <- max(heatmap_data_long\$row_group)  # Set row group to the current last group

      # Add the new padding rows to heatmap_data_long
      heatmap_data_long <- bind_rows_base(heatmap_data_long, padding_data)
    }

    return(heatmap_data_long)
  }

  # Load the CSV data
  heatmap_data <- read.csv(input_csv_path)

  # Check if the necessary column exists in the data
  if (!"total_counts_per_cov" %in% colnames(heatmap_data)) {
    stop("The column 'total_counts_per_cov' is not found in the data.")
  }

  # Create heatmap_data_long by selecting necessary columns
  heatmap_data_long <- heatmap_data[, c("mut_aa", "position", "total_counts_per_cov", "wt_aa")]  # Use 'total_counts_per_cov'

  # Find the smallest non-NA value in total_counts_per_cov
  min_non_na_value <- min(heatmap_data_long\$total_counts_per_cov, na.rm = TRUE)

  # Group positions by rows (75 positions per row) and calculate row_group
  heatmap_data_long\$row_group <- ((heatmap_data_long\$position - 1) %/% 75) + 1  # Grouping positions into rows

  # Apply padding to add missing positions at the end of the last row, using the calculated min value
  heatmap_data_long <- pad_heatmap_data_long(heatmap_data_long, min_non_na_value)

  # Convert positions to numeric, sort them, and create wt_aa_pos for the plot
  heatmap_data_long\$position <- as.numeric(heatmap_data_long\$position)  # Ensure position is numeric
  heatmap_data_long <- arrange_base(heatmap_data_long, heatmap_data_long\$position)  # Sort by position
  heatmap_data_long\$wt_aa_pos <- factor(paste0(heatmap_data_long\$wt_aa, heatmap_data_long\$position),
                                       levels = unique(paste0(heatmap_data_long\$wt_aa, heatmap_data_long\$position)))  # Create sorted factor levels for wt_aa_pos

  # Add a column to identify synonymous mutations (where mut_aa == wt_aa)
  heatmap_data_long\$synonymous <- heatmap_data_long\$mut_aa == heatmap_data_long\$wt_aa

  # Definiere die korrekte Reihenfolge der Aminosäuren
  amino_acid_order <- rev(c("G", "A", "V", "L", "M", "I", "F",
                         "Y", "W", "K", "R", "H", "D", "E",
                         "S", "T", "C", "N", "Q", "P", "*"))

  heatmap_data_long\$mut_aa <- factor(heatmap_data_long\$mut_aa, levels = amino_acid_order)

  # Bearbeite heatmap_data_long und erstelle syn_positions gleichzeitig
  syn_positions <- heatmap_data_long
  syn_positions\$mut_aa <- factor(syn_positions\$mut_aa, levels = amino_acid_order)
  # Berechne die x-Koordinate, die pro Gruppe immer von 1 bis 75 verläuft
  syn_positions\$x <- as.numeric(factor(syn_positions\$wt_aa_pos, levels = unique(syn_positions\$wt_aa_pos))) - ((syn_positions\$row_group - 1) * 75)
  syn_positions\$y <- as.numeric(factor(syn_positions\$mut_aa, levels = amino_acid_order))
  syn_positions <- filter_base(syn_positions, syn_positions\$synonymous == TRUE)

  # Calculate the number of row groups and adjust plot height dynamically
  num_row_groups <- max(heatmap_data_long\$row_group)
  plot_height <- num_row_groups * 4

  # Set the limits for the color scale, ignoring NA (negative values are now NA)
  min_count <- min(heatmap_data_long\$total_counts_per_cov, na.rm = TRUE)
  max_count <- max(heatmap_data_long\$total_counts_per_cov, na.rm = TRUE)
  max_position <- max(heatmap_data\$position)

  # Create the heatmap plot with explicit handling for positions > max_position
  heatmap_plot <- ggplot(heatmap_data_long, aes(x = wt_aa_pos, y = mut_aa, fill = total_counts_per_cov)) +
    scale_fill_gradientn(colours = c(alpha("blue", 0), "blue"), na.value = "grey35", trans = "log",  # Apply log transformation to the scale
                         limits = c(min_count, max_count),
                         breaks = scales::trans_breaks("log10", function(x) 10^x),  # Logarithmic scale breaks
                         labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_x_discrete(labels = function(x) {
      numeric_pos <- as.numeric(gsub("[^0-9]", "", x))
      ifelse(numeric_pos > max_position, " ", x)
    }) +
    geom_tile() +

    # Add diagonal lines for synonymous mutations using geom_segment
    geom_segment(data = syn_positions[syn_positions\$position <= max_position, ],
                 aes(x = x - 0.485, xend = x + 0.485,
                     y = y - 0.485, yend = y + 0.485, color = synonymous),
                 size = 0.2) +

    # Manuelle Farbskala für die diagonalen Linien
    scale_color_manual(values = c("TRUE" = "grey10"), labels = c("TRUE" = "")) +

    theme_minimal() +
    labs(title = "Heatmap of Counts per Coverage for Mutations", x = "Wild-type Amino Acid", y = "Mutant Amino Acid", fill = "Counts per \\n Coverage", color = "Synonymous Mutation") +
    theme(plot.title = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),  # Larger y-axis labels
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.title = element_text(size = 14),  # Larger legend title
          legend.text = element_text(size = 12),  # Larger legend text
          panel.spacing = unit(0.1, "lines"),  # Adjust panel spacing
          strip.text = element_blank(),  # Remove row group labels (facet numbers)
          strip.background = element_blank(),
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank()) +  # Remove minor grid lines
    facet_wrap(~ row_group, scales = "free_x", ncol = 1) +  # Group by 75 positions per row
    theme(panel.spacing = unit(0.2, "lines"))

  heatmap_plot <- heatmap_plot +
    geom_point(data = heatmap_data_long, aes(size = ""), colour = "black", alpha = 0)  # Invisible points for legend
  heatmap_plot <- heatmap_plot +
    guides(size = guide_legend(paste("Dropout (Counts <", threshold, ")"), override.aes = list(shape = 15, size = 8, colour = "grey35", alpha = 1)))  # Define Legend for Dropouts

  # Save the heatmap plot
  if (img_format == "pdf") {
    ggsave(output_pdf_path, plot = heatmap_plot, width = 16, height = plot_height, dpi = 150, device = cairo_pdf)
  } else {
    ggsave(output_pdf_path, plot = heatmap_plot, width = 16, height = plot_height, dpi = 150)
  }

  if (file.exists(output_pdf_path)) {
    print("Heatmap image successfully created!")
  } else {
    print("Error: Heatmap image was not created.")
  }
}

####################################################################################################
# 3. Sliding-window counts and counts per coverage along the ORF (rolling_counts.pdf, rolling_counts_per_cov.pdf)
####################################################################################################
# input: prefiltered gatk path, aa seq path, window_size (sliding window), output_path_folder
# output: two lineplots showing counts & counts per coverage divided in types of variants (single-/double-/triple-base exchange)

position_biases_counts <- function(prefiltered_gatk_path, aa_seq_path, window_size = 10) {

  # Load and process the data
  prefiltered_gatk <- read.table(prefiltered_gatk_path, sep = ",", fill = NA, header = TRUE)
  prefiltered_gatk\$pos <- as.numeric(sub("(\\\\D)(\\\\d+)(\\\\D)", "\\\\2", prefiltered_gatk\$pos_mut))
  unique_pos <- unique(as.numeric(prefiltered_gatk\$pos))
  aa_seq <- readLines(aa_seq_path, warn = FALSE)
  aa_seq_length <- nchar(aa_seq)
  aa_positions <- seq(nchar(aa_seq))

  means_counts_single <- rep(NA, nchar(aa_seq))
  means_counts_double <- rep(NA, nchar(aa_seq))
  means_counts_triple <- rep(NA, nchar(aa_seq))
  means_counts_per_cov_single <- rep(NA, nchar(aa_seq))
  means_counts_per_cov_double <- rep(NA, nchar(aa_seq))
  means_counts_per_cov_triple <- rep(NA, nchar(aa_seq))


  # The rows of each position, in table order (one pass instead of one scan of the table per position)
  rows_by_position <- split(seq_len(nrow(prefiltered_gatk)),
                            factor(match(prefiltered_gatk\$pos, aa_positions), levels = seq_along(aa_positions)))

  # Loop through each position in the amino acid sequence
  for (i in 1:(nchar(aa_seq))) {

    # Filter the data for the current position in aa_positions
    window_data <- prefiltered_gatk[rows_by_position[[i]], , drop = FALSE]

    # Calculate mean for Single mutations (where varying_bases == 1)
    window_data_single <- filter_base(window_data, window_data\$varying_bases == 1)
    means_counts_single[i] <- mean(window_data_single\$counts, na.rm = FALSE)
    means_counts_per_cov_single[i] <- mean(window_data_single\$counts_per_cov, na.rm = FALSE)

    # Calculate mean for Double mutations (where varying_bases == 2)
    window_data_double <- filter_base(window_data, window_data\$varying_bases == 2)
    means_counts_double[i] <- mean(window_data_double\$counts, na.rm = FALSE)
    means_counts_per_cov_double[i] <- mean(window_data_double\$counts_per_cov, na.rm = FALSE)

    # Calculate mean for Triple mutations (where varying_bases == 3)
    window_data_triple <- filter_base(window_data, window_data\$varying_bases == 3)
    means_counts_triple[i] <- mean(window_data_triple\$counts, na.rm = FALSE)
    means_counts_per_cov_triple[i] <- mean(window_data_triple\$counts_per_cov, na.rm = FALSE)
  }

  pos_bias_df <- data.frame(pos = seq(nchar(aa_seq)))


  pos_bias_df\$rolling_mean_counts_single <- roll_extend(means_counts_single, width = window_size,
                                               FUN = function(x) (mean(x, na.rm = TRUE) + 0.00001))
  pos_bias_df\$rolling_mean_counts_double <- roll_extend(means_counts_double, width = window_size,
                                                      FUN = function(x) (mean(x, na.rm = TRUE) + 0.00001))
  pos_bias_df\$rolling_mean_counts_triple <- roll_extend(means_counts_triple, width = window_size,
                                                      FUN = function(x) (mean(x, na.rm = TRUE) + 0.00001))

  pos_bias_df\$rolling_mean_counts_per_cov_single <- roll_extend(means_counts_per_cov_single, width = window_size,
                                                       FUN = function(x) (mean(x, na.rm = TRUE) + 0.00001))
  pos_bias_df\$rolling_mean_counts_per_cov_double <- roll_extend(means_counts_per_cov_double, width = window_size,
                                                              FUN = function(x) (mean(x, na.rm = TRUE) + 0.00001))
  pos_bias_df\$rolling_mean_counts_per_cov_triple <- roll_extend(means_counts_per_cov_triple, width = window_size,
                                                              FUN = function(x) (mean(x, na.rm = TRUE) + 0.00001))


  # Replace NAs with 0 for rolling means and SE
  pos_bias_df\$rolling_mean_counts_single[is.na(pos_bias_df\$rolling_mean_counts_single)] <- 0.00001
  pos_bias_df\$rolling_mean_counts_double[is.na(pos_bias_df\$rolling_mean_counts_double)] <- 0.00001
  pos_bias_df\$rolling_mean_counts_triple[is.na(pos_bias_df\$rolling_mean_counts_triple)] <- 0.00001
  pos_bias_df\$rolling_mean_counts_per_cov_single[is.na(pos_bias_df\$rolling_mean_counts_per_cov_single)] <- 0.000001
  pos_bias_df\$rolling_mean_counts_per_cov_double[is.na(pos_bias_df\$rolling_mean_counts_per_cov_double)] <- 0.000001
  pos_bias_df\$rolling_mean_counts_per_cov_triple[is.na(pos_bias_df\$rolling_mean_counts_per_cov_triple)] <- 0.000001



  plots_theme <- list(

    # Add the minimal theme to the list
    theme_minimal(),

    # Customize legend title and appearance
    theme(legend.title = element_text(size = 10, face = "bold"),
          legend.position = "right"),

    # Customize guides for the legend elements
    guides(
      color = guide_legend(title = "Mutation Type", order = 1),  # Title for line color
      fill = guide_legend(title = "Standard Deviation", order = 2),  # Title for ribbon fill
      linetype = guide_legend(title = "Required Coverage", order = 3)  # Title for line type
    )
  )


   # Placeholder name for the linetype
  linetype_placeholder <- "Required coverage"

  rolling_counts_plot <- ggplot(pos_bias_df, aes(x = pos)) +

    # Add line for rolling mean coverage
    geom_line(aes(y = rolling_mean_counts_single, color = "One Varying Base")) +
    geom_line(aes(y = rolling_mean_counts_double, color = "Two Varying Bases")) +
    geom_line(aes(y = rolling_mean_counts_triple, color = "Three Varying Bases")) +

    # Individual axis labels for this plot
    xlab("Amino Acid Position") +
    ylab("Counts") +

    # Manually set color and fill labels
    scale_color_manual(name = "Mutation Type",
                       values = c("One Varying Base" = "chocolate", "Two Varying Bases" = "darkolivegreen3", "Three Varying Bases" = "deepskyblue1"),
                       limits = c("One Varying Base", "Two Varying Bases", "Three Varying Bases")) + # Color for line

    # Apply the saved theme and design elements at the end
    plots_theme +

    scale_y_continuous(trans = 'log10', labels = scales::comma)




  rolling_counts_per_cov_plot <- ggplot(pos_bias_df, aes(x = pos)) +

    # Add line for rolling mean coverage
    geom_line(aes(y = rolling_mean_counts_per_cov_single, color = "One Varying Base")) +
    geom_line(aes(y = rolling_mean_counts_per_cov_double, color = "Two Varying Bases")) +
    geom_line(aes(y = rolling_mean_counts_per_cov_triple, color = "Three Varying Bases")) +

    # Individual axis labels for this plot
    xlab("Amino Acid Position") +
    ylab("Counts") +

    # Manually set color and fill labels
    scale_color_manual(name = "Mutation Type",
                       values = c("One Varying Base" = "chocolate", "Two Varying Bases" = "darkolivegreen3", "Three Varying Bases" = "deepskyblue1"),
                       limits = c("One Varying Base", "Two Varying Bases", "Three Varying Bases")) + # Color and legend order for lines

    # Apply the saved theme and design elements at the end
    plots_theme +

    scale_y_continuous(trans = 'log10', labels = scales::comma)


  ggsave(filename = "rolling_counts.pdf", plot = rolling_counts_plot, device = "pdf", width = 10, height = 6)
  ggsave(filename = "rolling_counts_per_cov.pdf", plot = rolling_counts_per_cov_plot, device = "pdf", width = 10, height = 6)
}

####################################################################################################
# 4. Sliding-window coverage along the ORF (rolling_coverage.pdf)
####################################################################################################
# input: prefiltered gatk path, aa seq path, window_size (sliding window), output_path_folder, aimed counts per aa variant
# output: lineplot showing coverage over positions and dotted line: Assuming all potential 21 variants are equally distributed, this is the coverage one need to get the aimed counts per variant.

position_biases_cov <- function(prefiltered_gatk_path, aa_seq_path, window_size = 10, output_file_path, targeted_counts_per_aa_variant = 100) {

  # Load and process the data
  prefiltered_gatk <- read.table(prefiltered_gatk_path, sep = ",", fill = NA, header = TRUE)
  prefiltered_gatk\$pos <- as.numeric(sub("(\\\\D)(\\\\d+)(\\\\D)", "\\\\2", prefiltered_gatk\$pos_mut))
  unique_pos <- unique(as.numeric(prefiltered_gatk\$pos))
  aa_seq <- readLines(aa_seq_path, warn = FALSE)
  aa_seq_length <- nchar(aa_seq)
  aa_positions <- seq(nchar(aa_seq))

  means_cov <- rep(NA, nchar(aa_seq))

  # Calculate means for cov (should be the same over all variants in one position)
  rows_by_position <- split(seq_len(nrow(prefiltered_gatk)),
                            factor(match(prefiltered_gatk\$pos, aa_positions), levels = seq_along(aa_positions)))
  for (i in 1:(nchar(aa_seq))) {
    window_data <- prefiltered_gatk[rows_by_position[[i]], , drop = FALSE]
    means_cov[i] <- mean(window_data\$cov, na.rm = FALSE)
  }

  pos_bias_df <- data.frame(pos = seq(nchar(aa_seq)))

  # Log-transform the rolling means to avoid issues with zeros (log(y + 0.001))
  pos_bias_df\$rolling_mean_cov <- roll_extend(means_cov, width = window_size,
                                            FUN = function(x) (mean(x, na.rm = TRUE)))

  # Generate the plot
  plots_theme <- list(

    # Customize legend appearance (leave titles blank)
    guides(
      color = guide_legend(title = NULL, order = 1),  # Title for line color
      fill = guide_legend(title = NULL, order = 2),  # Title for ribbon fill
      linetype = guide_legend(title = NULL, order = 3)  # Title for line type
    ),

    # Apply minimal theme
    theme_minimal()
  )


  linetype_placeholder <- "Required coverage"

  rolling_cov_plot <- ggplot(pos_bias_df, aes(x = pos, y = rolling_mean_cov)) +

    # Add line for rolling mean coverage
    geom_line(aes(color = "Coverage")) +

    # Add a horizontal line with a mapped linetype so it appears in the legend
    geom_hline(aes(yintercept = ((21 * nchar(aa_seq) * targeted_counts_per_aa_variant)),
                   linetype = linetype_placeholder),
               color = "black", linewidth = 0.3) +

    # Individual axis labels for this plot
    xlab("Amino Acid Position") +
    ylab("Coverage") +

    # Manually set color and fill labels
    scale_color_manual(values = c("Coverage" = "black")) +  # Color for line

    # Set the linetype with a label that uses the targeted threshold dynamically
    scale_linetype_manual(values = c("Required coverage" = "dotted"),
                          labels = paste("Required coverage \\n for", as.character(targeted_counts_per_aa_variant), "counts per variant \\n if equally present")) +

    # Apply the saved theme and design elements at the end
    plots_theme +

    # Set y-axis to log scale and apply comma formatting
    scale_y_continuous(trans = 'log10', labels = scales::comma)


  # Return the plot
  ggsave(filename = output_file_path, plot = rolling_cov_plot, device = "pdf", width = 10, height = 6)
}

####################################################################################################
# 5. LogDiff: sorted counts per coverage of all programmed variants (logdiff_plot.pdf, logdiff_varying_bases.pdf)
####################################################################################################

# input: completed library table (all programmed variants; counts_per_cov 0 for variants with no
#        counts, observed = FALSE for those never observed)
# output: two LogDiff plots (counts per coverage, sorted, log scale): a line plot, and a dot plot
#         coloured by the number of changed bases per codon
#
# Variants with no counts cannot be drawn on a log scale. They are placed at a floor derived from
# the data - half the lowest non-zero counts per coverage in the library, i.e. below every observed
# variant - and shown in grey. LogDiff is log10 of the 90th over the 10th percentile of all
# programmed variants. When more than 10 % of the library has no counts, the 10th percentile would
# be the floor, an artefact of the plot, so LogDiff is not reported; the share of variants without
# counts, printed on every plot, is then the meaningful measure.

logdiff_plot <- function(completed_path) {

  completed <- read.table(completed_path, sep = ",", fill = NA, header = TRUE)
  cpc <- completed\$counts_per_cov
  cpc[is.na(cpc)] <- 0
  if (!"observed" %in% names(completed)) {
    completed\$observed <- cpc > 0
  }

  n_total    <- length(cpc)
  n_never    <- sum(!completed\$observed)
  n_zeroed   <- sum(completed\$observed & cpc == 0)    # observed, but corrected to 0
  n_no_count <- n_never + n_zeroed

  nonzero <- cpc[cpc > 0]
  if (length(nonzero) == 0) {
    # Nothing to rank: write placeholder plots that say so, instead of failing the run
    empty <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No programmed variant has counts in this library.") +
      theme_void()
    ggsave(filename = "logdiff_plot.pdf", plot = empty, device = "pdf", width = 10, height = 6)
    ggsave(filename = "logdiff_varying_bases.pdf", plot = empty, device = "pdf", width = 10, height = 6)
    return(invisible(NULL))
  }
  floor_value <- min(nonzero) / 2

  sorted_counts <- completed
  sorted_counts\$counts_per_cov_plot <- ifelse(cpc > 0, cpc, floor_value)
  sorted_counts\$category <- ifelse(cpc > 0, as.character(sorted_counts\$varying_bases), "none")
  sorted_counts <- arrange_base(sorted_counts, sorted_counts\$counts_per_cov_plot)
  sorted_counts\$ids <- seq_len(nrow(sorted_counts))

  Q3 <- quantile(sorted_counts\$counts_per_cov_plot, 0.90, na.rm = TRUE)
  if (n_no_count / n_total > 0.10) {
    Q1 <- NA
    logdiff_text <- "LogDiff: not defined (more than 10 % of variants without counts)"
  } else {
    Q1 <- quantile(sorted_counts\$counts_per_cov_plot, 0.10, na.rm = TRUE)
    logdiff_text <- paste0("LogDiff = ", round(log10(Q3) - log10(Q1), 2))
  }
  percentile_lines <- c(Q1, Q3)[!is.na(c(Q1, Q3))]

  label <- paste0(
    logdiff_text, "\\n",
    "No counts: ", n_no_count, " of ", n_total, " variants (", round(100 * n_no_count / n_total, 1), " %)",
    if (n_zeroed > 0) paste0(", ", n_zeroed, " of them corrected to 0") else ""
  )
  caption <- paste0(
    "Variants without counts are drawn at half the lowest observed counts per coverage (",
    format(signif(floor_value, 3), scientific = TRUE), ", dashed line)."
  )
  y_top <- max(sorted_counts\$counts_per_cov_plot)

  line_plot <- ggplot(sorted_counts, aes(x = ids, y = counts_per_cov_plot)) +
    geom_line(color = "black") +
    geom_hline(yintercept = floor_value, linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = percentile_lines, linetype = "dotted", color = "black") +
    annotate("text", x = 0, y = y_top, label = label, hjust = 0, vjust = 1, size = 4.5, color = "black") +
    scale_y_continuous(trans = "log10") +
    labs(x = "Variants", y = "Counts per Coverage", caption = caption) +
    theme_minimal()
  ggsave(filename = "logdiff_plot.pdf", plot = line_plot, device = "pdf", width = 10, height = 6)

  colored_plot <- ggplot(sorted_counts, aes(x = ids, y = counts_per_cov_plot, color = category)) +
    geom_hline(yintercept = floor_value, linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = percentile_lines, linetype = "dotted", color = "black") +
    geom_point(size = 0.9) +
    annotate("text", x = 1, y = y_top, label = label, hjust = 0, vjust = 1, size = 4.5, color = "black") +
    scale_y_continuous(trans = "log10") +
    scale_color_manual(values = c("1" = "chocolate", "2" = "darkolivegreen3", "3" = "deepskyblue1", "none" = "grey70"),
                       breaks = c("1", "2", "3", "none"),
                       labels = c("1", "2", "3", "no counts"),
                       name = "Varying \\n Bases") +
    labs(x = "Variants", y = "Counts per Coverage", caption = caption) +
    theme_minimal()
  ggsave(filename = "logdiff_varying_bases.pdf", plot = colored_plot, device = "pdf", width = 10, height = 6)
}

####################################################################################################
# 6. Sequencing-depth rarefaction, with --run_seqdepth (SeqDepth.pdf, seqdepth_curve.csv)
####################################################################################################
# Sequencing-depth rarefaction: what fraction of the possible library variants stay detectable
# (>= `threshold` reads) as the sequencing depth is subsampled from full down to zero.
#
# This is a closed-form rarefaction, not a stochastic simulation. Subsampling reads without
# replacement makes each variant's surviving read count hypergeometric, so the expected number of
# variants still at or above the threshold at depth n is a sum of hypergeometric tail probabilities:
#
#   E[detected | n reads] = sum_i  P(Hypergeom(N, c_i, n) >= threshold)
#                         = sum_i  phyper(threshold - 1, c_i, N - c_i, n, lower.tail = FALSE)
#
# where N is the total read count and c_i the reads on variant i. It is exact (no Monte-Carlo noise),
# needs no random sampling, and is orders of magnitude faster than draw-one-read-at-a-time loops - so
# the step count can be fine without the runtime blowing up. Being exact at every depth, a finer grid
# only resolves the curve better; it cannot change it.
#
# Depth grid: the step size grows by a constant factor from 0.01% at the bottom to 5% at the top
# (0.95 -> 1). Each step is then a fixed fraction (~5%) of the depth it starts from, so the steep
# low-depth rise is resolved as finely as the flat top, whatever depth a given library saturates at.
# With these defaults that is 123 steps.

# 0, then the cumulative sum of K geometrically growing steps that ends at exactly 1 with a final step
# of `last_step`. K comes from the closed-form ratio (1 - first_step) / (1 - last_step), rounded up to
# a whole number of steps; the ratio is then re-solved for that K so both ends stay exact.
depth_grid <- function(first_step = 1e-4, last_step = 0.05) {
  K <- ceiling(1 + log(last_step / first_step) / log((1 - first_step) / (1 - last_step)))
  r <- uniroot(function(r) last_step * (1 - r^-K) / (1 - 1 / r) - 1,
               c(1 + 1e-9, 2), tol = 1e-12)\$root
  c(0, cumsum(last_step * r^-((K - 1):0)))
}

SeqDepth_simulation_plot <- function(prefiltered_gatk_path, possible_mutations_path,
                                     output_file_path, curve_csv_path,
                                     first_step = 1e-4, last_step = 0.05,
                                     threshold = 3) {

  data <- read.csv(prefiltered_gatk_path)
  counts <- as.numeric(data\$counts)
  counts <- counts[!is.na(counts) & counts > 0]           # only observed variants carry reads
  possible_mutations <- read.csv(possible_mutations_path)
  total_possible <- nrow(possible_mutations)               # denominator: the whole designed library

  N <- sum(counts)
  if (N <= 0 || total_possible <= 0) {
    stop("SeqDepth: no counts or no possible mutations to rarefy")
  }

  # Depth fractions from depth_grid(), as integer read totals. unique() drops the collisions that
  # rounding produces at the bottom of very small libraries.
  depths <- unique(round(depth_grid(first_step, last_step) * N))

  # Expected detected variants at each depth. phyper is vectorised over the per-variant counts, so each
  # depth is one vectorised call; `lower.tail = FALSE` gives P(X > threshold-1) = P(X >= threshold).
  expected_detected <- vapply(depths, function(n) {
    if (n <= 0) return(0)
    sum(phyper(threshold - 1, counts, N - counts, n, lower.tail = FALSE))
  }, numeric(1))

  curve <- data.frame(
    depth_fold   = round(depths / N, 6),                       # x: fraction of full sequencing depth
    variants_pct = round(expected_detected / total_possible * 100, 4)  # y: % of the possible library
  )
  write.csv(curve, curve_csv_path, row.names = FALSE)

  plot <- ggplot(curve, aes(x = depth_fold, y = variants_pct)) +
    geom_line(color = "black", linewidth = 0.4) +
    geom_hline(yintercept = 100, linetype = "dotted", color = "black") +
    scale_y_continuous(
      labels = scales::percent_format(scale = 1),
      limits = c(0, 100),
      breaks = seq(0, 100, by = 5)
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, length.out = 20),
      labels = scales::number_format(accuracy = 0.01)
    ) +
    labs(x = "Fraction of sequencing depth", y = "Variants detected (% of possible library)") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid.major = element_line(linewidth = 0.2, linetype = "solid", color = "grey80"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(output_file_path, plot = plot, device = "pdf", width = 8, height = 6)
}

#####
# run functions
#####
# A plot that cannot be drawn for this library - typically one with too few variants, or variants at too
# few positions for the sliding window - is skipped: a warning in the task log, and a page saying why in
# each PDF it would have written, so that the outputs, the run report and a resumed run stay complete.
qc_step <- function(name, pdfs, expr) {
  tryCatch(expr, error = function(e) {
    msg <- paste0(name, " not drawn for this library: ", conditionMessage(e))
    message("WARNING: ", msg)
    while (!is.null(grDevices::dev.list())) grDevices::dev.off()
    for (f in pdfs) {
      grDevices::pdf(f, width = 10, height = 6)
      graphics::plot.new()
      graphics::text(0.5, 0.5, paste(strwrap(msg, 90), collapse = "\\n"))
      grDevices::dev.off()
    }
    invisible(NULL)
  })
}

qc_step("Counts heatmap", "counts_heatmap.pdf", counts_heatmap(
  input_csv_path = "$variantCounts_for_heatmaps",
  threshold = 1,
  output_pdf_path = "counts_heatmap.pdf",
  img_format = "pdf"
))

qc_step("Counts-per-coverage heatmap", "counts_per_cov_heatmap.pdf", counts_per_cov_heatmap(
  input_csv_path = "$variantCounts_for_heatmaps",
  threshold = 1,
  output_pdf_path = "counts_per_cov_heatmap.pdf",
  img_format = "pdf"
))

qc_step("Sliding-window counts", c("rolling_counts.pdf", "rolling_counts_per_cov.pdf"), position_biases_counts(
  prefiltered_gatk_path = "$variantCounts_filtered_by_library",
  aa_seq_path = "$aa_seq",
  window_size = $sliding_window_size
))

qc_step("Sliding-window coverage", "rolling_coverage.pdf", position_biases_cov(
  prefiltered_gatk_path = "$variantCounts_filtered_by_library",
  aa_seq_path = "$aa_seq",
  window_size = $sliding_window_size,
  output_file_path = "rolling_coverage.pdf",
  targeted_counts_per_aa_variant = $aimed_cov
))

qc_step("LogDiff", c("logdiff_plot.pdf", "logdiff_varying_bases.pdf"), logdiff_plot(
  completed_path = "$library_completed_variantCounts"
))

if (as.logical("$run_seqdepth")) {
  qc_step("Sequencing-depth rarefaction", "SeqDepth.pdf", SeqDepth_simulation_plot(
    prefiltered_gatk_path   = "$variantCounts_filtered_by_library",
    possible_mutations_path = "$possible_mutations",
    output_file_path        = "SeqDepth.pdf",
    curve_csv_path          = "seqdepth_curve.csv",
    first_step              = 1e-4,
    last_step               = 0.05,
    threshold               = 1
  ))
}

####
# create versions.yml
####
r_version <- strsplit(version[['version.string']], ' ')[[1]][3]
ggplot2_version <- as.character(packageVersion("ggplot2"))
scales_version <- as.character(packageVersion("scales"))

if (is.null(r_version)) r_version <- "unknown"
if (length(ggplot2_version) == 0) ggplot2_version <- "unknown"
if (length(scales_version) == 0) scales_version <- "unknown"

f <- file("versions.yml", "w")
writeLines(
  c(
    '"${task.process}":',
    paste('    r-base:', r_version),
    paste('    r-ggplot2:', ggplot2_version),
    paste('    r-scales:', scales_version)
  ),
  f
)
close(f)
