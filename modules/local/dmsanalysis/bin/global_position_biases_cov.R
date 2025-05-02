# input: prefiltered gatk path, aa seq path, window_size (sliding window), output_path_folder, aimed counts per aa variant
# output: lineplot showing coverage over positions and dotted line: Assuming all potential 21 variants are equally distributed, this is the coverage one need to get the aimed counts per variant.

library(zoo) # sliding window
library(dplyr)
library(ggplot2)
library(scales)
position_biases <- function(prefiltered_gatk_path, aa_seq_path, window_size = 10, output_file_path, targeted_counts_per_aa_variant = 100) {

  # Load and process the data
  prefiltered_gatk <- read.table(prefiltered_gatk_path, sep = ",", fill = NA, header = TRUE)
  prefiltered_gatk$pos <- as.numeric(sub("(\\D)(\\d+)(\\D)", "\\2", prefiltered_gatk$pos_mut))
  unique_pos <- unique(as.numeric(prefiltered_gatk$pos))
  aa_seq <- readLines(aa_seq_path, warn = FALSE)
  aa_seq_length <- nchar(aa_seq)
  aa_positions <- seq(nchar(aa_seq))

  means_cov <- rep(NA, nchar(aa_seq))

  # Calculate means for cov (should be the same over all variants in one position)
  for (i in 1:(nchar(aa_seq))) {
    window_data <- prefiltered_gatk %>% filter(prefiltered_gatk$pos %in% aa_positions[i])
    means_cov[i] <- mean(window_data$cov, na.rm = FALSE)
  }

  pos_bias_df <- data.frame(pos = seq(nchar(aa_seq)))

  # Log-transform the rolling means to avoid issues with zeros (log(y + 0.001))
  pos_bias_df$rolling_mean_cov <- rollapply(means_cov, width = window_size,
                                            FUN = function(x) (mean(x, na.rm = TRUE)), fill = "extend")

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
                          labels = paste("Required coverage \n for", as.character(targeted_counts_per_aa_variant), "counts per variant \n if equally present")) +

    # Apply the saved theme and design elements at the end
    plots_theme +

    # Set y-axis to log scale and apply comma formatting
    scale_y_continuous(trans = 'log10', labels = scales::comma)


  # Return the plot
  ggsave(filename = output_file_path, plot = rolling_cov_plot, device = "pdf", width = 10, height = 6)
}

# Example call to the function
#position_biases("/Users/benjaminwehnert/CRG/DMS_QC/testing_data/gatk_filtered_by_codon_library.csv", "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/aa_seq.txt", window_size = 18, output_path_folder = "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/testing_outputs/rolling_coverage.pdf", targeted_counts_per_aa_variant = 15)




