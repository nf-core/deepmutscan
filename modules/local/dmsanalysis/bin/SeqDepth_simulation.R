# input: prefiltered (by codon library) gatk path, possible mutations path, output_folder, reduction_fraction (steps in % to reduce counts), threshold to count variant as present in dataset
# output: pdf showing the plot
# limitation: takes quite a lot time: depends on number of counts in prefiltered gatk (4 min on M1 MacBook for 340,000 counts in total) -> reduction_fraction only has a low impact -> need to find more efficient random-sampling algorithm

library(dplyr)
library(ggplot2)

SeqDepth_simulation_plot <- function(prefiltered_gatk_path, possible_mutations_path, output_file_path, reduction_fraction = 0.01, threshold = 3) {

  # Read data from the specified CSV file
  data <- read.csv(prefiltered_gatk_path)
  data <- data %>% mutate(counts = as.numeric(counts)) # Ensure counts are numeric
  original_counts <- data$counts # Store the original counts for weight calculations
  possible_mutations <- read.csv(possible_mutations_path)

  # Initialize variables
  total_counts <- sum(data$counts)
  reduction_per_step <- floor(total_counts * reduction_fraction) # Round down to the nearest integer
  results <- data.frame(remaining_counts = numeric(), remaining_variants = numeric())

  # Track the initial state
  remaining_variants <- sum(data$counts >= threshold)
  remaining_counts <- total_counts
  results <- rbind(results, data.frame(remaining_counts = remaining_counts, remaining_variants = remaining_variants))

  # Get a list of indices for non-zero variants
  non_zero_indices <- which(data$counts > 0)

  # Calculate initial weights based on the original counts for non-zero variants
  weights <- as.numeric(original_counts[non_zero_indices])
  weights <- weights / sum(weights) # Normalize weights

  # Loop until all counts are zero
  while (remaining_counts > 0) {
    # If no more non-zero variants are available, break the loop
    if (length(non_zero_indices) == 0) {
      break
    }

    # Randomly reduce counts by the specified amount using weighted sampling
    for (i in 1:reduction_per_step) {
      # Randomly choose an index from the non-zero variants based on weights
      selected_idx <- sample(length(non_zero_indices), 1, prob = weights)
      index <- non_zero_indices[selected_idx]

      # Reduce the count by 1
      data$counts[index] <- data$counts[index] - 1

      # If the count reaches zero, remove the index from the list of non-zero variants
      if (data$counts[index] == 0) {
        non_zero_indices <- non_zero_indices[-selected_idx]
        weights <- weights[-selected_idx] # Remove the corresponding weight
      }
    }

    # Update remaining counts and remaining variants after the reduction
    remaining_counts <- sum(data$counts)
    remaining_variants <- sum(data$counts >= threshold & data$counts > 0)

    # Store the results for this step
    results <- rbind(results, data.frame(remaining_counts = remaining_counts, remaining_variants = remaining_variants))

    # Adjust the reduction_per_step if the total remaining counts are less than the reduction amount
    if (remaining_counts < reduction_per_step) {
      reduction_per_step <- remaining_counts
    }

    # Recalculate weights using the original counts, but only for remaining non-zero variants
    if (length(non_zero_indices) > 0) {
      weights <- as.numeric(original_counts[non_zero_indices])
      weights <- weights / sum(weights) # Normalize weights
    }
  }

  # Transform results for plotting
  baseline_count <- max(results$remaining_counts)
  plot_data <- results %>%
    mutate(
      remaining_counts_fold = round(remaining_counts / baseline_count, 2),  # X-axis in fold-change from max, rounded to 2 decimals
      remaining_variants_percent = (remaining_variants / nrow(possible_mutations)) * 100  # Y-axis in percent
    )

  # Set plot limits
  x_max <- 1  # X-axis limit exactly at 1
  y_max <- 100  # Y-axis limit at 100%

  # Create the plot
  plot <- ggplot(plot_data, aes(x = remaining_counts_fold, y = remaining_variants_percent)) +
    geom_line(color = "black", size = 0.4) +  # Solid line for the data
    geom_hline(yintercept = 100, linetype = "dotted", color = "black") +  # Horizontal line at 100%

    # Main plot settings with fine grid lines
    scale_y_continuous(
      labels = scales::percent_format(scale = 1),
      limits = c(0, y_max),
      breaks = seq(0, y_max, by = 5)  # Y-axis grid lines every 5%
    ) +
    scale_x_continuous(
      limits = c(0, x_max),
      breaks = seq(0, x_max, length.out = 20),  # X-axis with 20 grid lines ending at 1
      labels = scales::number_format(accuracy = 0.01)  # Round labels to 2 decimal places
    ) +
    labs(
      x = "Fold-Change of Sequencing Depth",
      y = "Variants (% of Maximum)"
    ) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA),  # Add a box-like border
      panel.grid.major = element_line(size = 0.2, linetype = "solid", color = "grey80"),  # Fine grid lines
      panel.grid.minor = element_blank(),  # No minor grid lines for clarity
      axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate X-axis labels at 45 degrees
    )

  # Save the plot as a PDF in the specified output folder
  ggsave(output_file_path, plot = plot, device = "pdf", width = 8, height = 6)
}

#Example usage with input and output paths
# results_weighted <- coverage_simulation_plot(
#   prefiltered_gatk_path = "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/gatk_filtered_by_codon_library.csv",
#   possible_mutations_path = "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/possible_mutations.csv",
#   output_folder = "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/testing_outputs",
#   reduction_fraction = 0.01,
#   threshold = 3
# )


#SeqDepth


