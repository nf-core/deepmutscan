
# Maybe develop additional ideas to characterise variants below 10 percentile -> Do we find patterns that help the user to understand why certain variants are hard to count?

# input: completed_prefiltered_gatk path, output folder path
# output: two logdiff-plots (counts_per_cov): 1st lineplot 2nd dotplot showing type of variant (varying bases - 1/2/3)

library(dplyr)
library(ggplot2)
library(scales)

 logdiff_plot <- function(prefiltered_gatk_path) {

   # Load the data
   prefiltered_gatk <- read.table(prefiltered_gatk_path, sep = ",", fill = NA, header = TRUE)

   # Sort by counts_per_cov while keeping corresponding varying_bases
   sorted_counts <- prefiltered_gatk %>%
     arrange(counts_per_cov)  # Sort by counts_per_cov

   # Create a new column for the sorted index (1, 2, 3, ...)
   sorted_counts$ids <- 1:nrow(sorted_counts)

   # Calculate Q1 (10% quantile) and Q3 (90% quantile)
   Q1 <- quantile(sorted_counts$counts_per_cov, 0.10, na.rm = TRUE)
   Q3 <- quantile(sorted_counts$counts_per_cov, 0.90, na.rm = TRUE)

   # Calculate the LogDiff
   LogDiff <- log10(Q3) - log10(Q1)

   # Create the first plot with a line
   line_plot <- ggplot(sorted_counts, aes(x = ids, y = counts_per_cov)) +
     geom_line(color = "black") +

     # Set axis labels
     xlab("Variants") +
     ylab("Counts per Coverage") +

     # Apply logarithmic scale to the y-axis
     scale_y_continuous(trans = 'log10') +

     # Add horizontal dotted lines at Q1 and Q3
     geom_hline(yintercept = Q1, linetype = "dotted", color = "black") +
     geom_hline(yintercept = Q3, linetype = "dotted", color = "black") +

     # Add the LogDiff value to the top left corner
     annotate("text", x = 0, y = max(sorted_counts$counts_per_cov), label = paste("LogDiff =", round(LogDiff, 2)), hjust = 0, vjust = 1, size = 5, color = "black") +

     # Apply the minimal theme
     theme_minimal()

   # Save the line plot as a PDF
   ggsave(filename = "logdiff_plot.pdf", plot = line_plot, device = "pdf", width = 10, height = 6)

   # Create the second plot with colored dots based on varying_bases
   colored_plot <- ggplot(sorted_counts, aes(x = ids, y = counts_per_cov, color = as.factor(varying_bases))) +

     # Add horizontal dotted lines at Q1 and Q3
     geom_hline(yintercept = Q1, linetype = "dotted", color = "black") +
     geom_hline(yintercept = Q3, linetype = "dotted", color = "black") +

     geom_point(size = 0.9) +  # Use points instead of lines

     # Set axis labels
     xlab("Variants") +
     ylab("Counts per Coverage") +

     # Apply logarithmic scale to the y-axis
     scale_y_continuous(trans = 'log10') +

     # Add the LogDiff value to the top left corner
     annotate("text", x = 1, y = max(sorted_counts$counts_per_cov), label = paste("LogDiff =", round(LogDiff, 2)), hjust = 0, vjust = 1, size = 5, color = "black") +

     # Add color scale for varying_bases
     scale_color_manual(values = c("1" = "chocolate", "2" = "darkolivegreen3", "3" = "deepskyblue1"), name = "Varying \n Bases") +

     # Apply the minimal theme
     theme_minimal()

   # Save the colored plot as a PDF
   ggsave(filename = "logdiff_varying_bases.pdf", plot = colored_plot, device = "pdf", width = 10, height = 6)
 }

 # Example call of the function
#  logdiff_plot("/Users/benjaminwehnert/CRG/DMS_QC/testing_data/testing_outputs/completed_prefiltered_gatk.csv", "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/testing_outputs")





 # logdiff_plot_codon_mut <- function(prefiltered_gatk_path, output_folder_path) {
 #
 #   # Load the data
 #   prefiltered_gatk <- read.table(prefiltered_gatk_path, sep = ",", fill = NA, header = TRUE)
 #
 #   # Extract the codon number from the codon_mut column (e.g., '377' from '377:TCC>CTG')
 #   prefiltered_gatk <- prefiltered_gatk %>%
 #     mutate(position = as.numeric(gsub(":.*", "", codon_mut)))  # Extract the number before the colon
 #
 #   # Sort by counts_per_cov while keeping corresponding position
 #   sorted_counts <- prefiltered_gatk %>%
 #     arrange(counts_per_cov)  # Sort by counts_per_cov
 #
 #   # Create a new column for the sorted index (1, 2, 3, ...)
 #   sorted_counts$ids <- 1:nrow(sorted_counts)
 #
 #   # Calculate Q1 (10% quantile) and Q3 (90% quantile)
 #   Q1 <- quantile(sorted_counts$counts_per_cov, 0.10, na.rm = TRUE)
 #   Q3 <- quantile(sorted_counts$counts_per_cov, 0.90, na.rm = TRUE)
 #
 #   # Calculate the LogDiff
 #   LogDiff <- log10(Q3) - log10(Q1)
 #
 #   # Create the second plot with a larger gray area and preserved horizontal lines and LogDiff
 #   codon_mut_plot <- ggplot(sorted_counts, aes(x = ids, y = counts_per_cov, color = position)) +
 #
 #     # Add horizontal dotted lines in the background
 #     geom_hline(yintercept = Q1, linetype = "dotted", color = "black") +
 #     geom_hline(yintercept = Q3, linetype = "dotted", color = "black") +
 #
 #     # Add points with smaller size
 #     geom_point(size = 0.9) +
 #
 #     # Set axis labels
 #     xlab("Sorted Index") +
 #     ylab("Counts per Coverage") +
 #
 #     # Apply logarithmic scale to the y-axis
 #     scale_y_continuous(trans = 'log10', labels = scales::comma) +
 #
 #     # Add the LogDiff value to the top left corner
 #     annotate("text", x = 1, y = max(sorted_counts$counts_per_cov), label = paste("LogDiff =", round(LogDiff, 2)), hjust = 0, vjust = 1, size = 5, color = "black") +
 #
 #     # Add a color gradient with a larger gray area
 #     scale_color_gradient2(low = "blue", mid = "lightgray", high = "green", midpoint = median(sorted_counts$position, na.rm = TRUE), name = "Position") +
 #
 #     # Apply the minimal theme
 #     theme_minimal()
 #
 #   # Save the second plot as a PDF
 #   ggsave(filename = paste(output_folder_path, "/logdiff_position.pdf", sep = ""),
 #          plot = codon_mut_plot, device = "pdf", width = 10, height = 6)
 # }
 #
 # # Example call
 # logdiff_plot_codon_mut("/Users/benjaminwehnert/CRG/DMS_QC/testing_data/testing_outputs/completed_prefiltered_gatk.csv", "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/testing_outputs")
 #



