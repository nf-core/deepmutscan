# input: prefiltered gatk path, aa seq path, window_size (sliding window), output_path_folder
# output: two lineplots showing counts & counts per coverage divided in types of variants (single-/double-/triple-base exchange)

library(zoo) # sliding window
library(dplyr)
library(ggplot2)
library(scales)

position_biases <- function(prefiltered_gatk_path, aa_seq_path, window_size = 10) {

  # Load and process the data
  prefiltered_gatk <- read.table(prefiltered_gatk_path, sep = ",", fill = NA, header = TRUE)
  prefiltered_gatk$pos <- as.numeric(sub("(\\D)(\\d+)(\\D)", "\\2", prefiltered_gatk$pos_mut))
  unique_pos <- unique(as.numeric(prefiltered_gatk$pos))
  aa_seq <- readLines(aa_seq_path, warn = FALSE)
  aa_seq_length <- nchar(aa_seq)
  aa_positions <- seq(nchar(aa_seq))

  means_counts_single <- rep(NA, nchar(aa_seq))
  means_counts_double <- rep(NA, nchar(aa_seq))
  means_counts_triple <- rep(NA, nchar(aa_seq))
  means_counts_per_cov_single <- rep(NA, nchar(aa_seq))
  means_counts_per_cov_double <- rep(NA, nchar(aa_seq))
  means_counts_per_cov_triple <- rep(NA, nchar(aa_seq))


  # Loop through each position in the amino acid sequence
  for (i in 1:(nchar(aa_seq))) {

    # Filter the data for the current position in aa_positions
    window_data <- prefiltered_gatk %>% filter(prefiltered_gatk$pos %in% aa_positions[i])

    # Calculate mean for Single mutations (where varying_bases == 1)
    window_data_single <- window_data %>% filter(varying_bases == 1)
    means_counts_single[i] <- mean(window_data_single$counts, na.rm = FALSE)
    means_counts_per_cov_single[i] <- mean(window_data_single$counts_per_cov, na.rm = FALSE)

    # Calculate mean for Double mutations (where varying_bases == 2)
    window_data_double <- window_data %>% filter(varying_bases == 2)
    means_counts_double[i] <- mean(window_data_double$counts, na.rm = FALSE)
    means_counts_per_cov_double[i] <- mean(window_data_double$counts_per_cov, na.rm = FALSE)

    # Calculate mean for Triple mutations (where varying_bases == 3)
    window_data_triple <- window_data %>% filter(varying_bases == 3)
    means_counts_triple[i] <- mean(window_data_triple$counts, na.rm = FALSE)
    means_counts_per_cov_triple[i] <- mean(window_data_triple$counts_per_cov, na.rm = FALSE)
  }

  pos_bias_df <- data.frame(pos = seq(nchar(aa_seq)))


  pos_bias_df$rolling_mean_counts_single <- rollapply(means_counts_single, width = window_size,
                                               FUN = function(x) (mean(x, na.rm = TRUE) + 0.00001), fill = "extend")
  pos_bias_df$rolling_mean_counts_double <- rollapply(means_counts_double, width = window_size,
                                                      FUN = function(x) (mean(x, na.rm = TRUE) + 0.00001), fill = "extend")
  pos_bias_df$rolling_mean_counts_triple <- rollapply(means_counts_triple, width = window_size,
                                                      FUN = function(x) (mean(x, na.rm = TRUE) + 0.00001), fill = "extend")

  pos_bias_df$rolling_mean_counts_per_cov_single <- rollapply(means_counts_per_cov_single, width = window_size,
                                                       FUN = function(x) (mean(x, na.rm = TRUE) + 0.00001), fill = "extend")
  pos_bias_df$rolling_mean_counts_per_cov_double <- rollapply(means_counts_per_cov_double, width = window_size,
                                                              FUN = function(x) (mean(x, na.rm = TRUE) + 0.00001), fill = "extend")
  pos_bias_df$rolling_mean_counts_per_cov_triple <- rollapply(means_counts_per_cov_triple, width = window_size,
                                                              FUN = function(x) (mean(x, na.rm = TRUE) + 0.00001), fill = "extend")


  # Replace NAs with 0 for rolling means and SE
  pos_bias_df$rolling_mean_counts_single[is.na(pos_bias_df$rolling_mean_counts_single)] <- 0.00001
  pos_bias_df$rolling_mean_counts_double[is.na(pos_bias_df$rolling_mean_counts_double)] <- 0.00001
  pos_bias_df$rolling_mean_counts_triple[is.na(pos_bias_df$rolling_mean_counts_triple)] <- 0.00001
  pos_bias_df$rolling_mean_counts_per_cov_single[is.na(pos_bias_df$rolling_mean_counts_per_cov_single)] <- 0.000001
  pos_bias_df$rolling_mean_counts_per_cov_double[is.na(pos_bias_df$rolling_mean_counts_per_cov_double)] <- 0.000001
  pos_bias_df$rolling_mean_counts_per_cov_triple[is.na(pos_bias_df$rolling_mean_counts_per_cov_triple)] <- 0.000001



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

# Example call to the function
#position_biases("/Users/benjaminwehnert/CRG/DMS_QC/testing_data/gatk_filtered_by_codon_library.csv", "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/aa_seq.txt", window_size = 18, output_path_folder = "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/testing_outputs")














