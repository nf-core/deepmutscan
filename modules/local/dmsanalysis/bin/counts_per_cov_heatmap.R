# Input: prefiltered GATK data path, output_path, threshold (same as used for prepare_gatk_data_for_counts_per_cov_heatmap function)
# Output: counts_per_cov_heatmap.pdf

library(dplyr)
library(ggplot2)

counts_per_cov_heatmap <- function(input_csv_path, threshold = 3, output_pdf_path, img_format = "pdf") {

  # Inner function to add padding to the last row, adding 21 amino acids per position
  pad_heatmap_data_long <- function(heatmap_data_long, min_non_na_value, num_positions_per_row = 75) {
    all_amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                         "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "*")

    max_position <- max(heatmap_data_long$position)
    num_missing_positions <- num_positions_per_row - (max_position %% num_positions_per_row)

    if (num_missing_positions < num_positions_per_row) {
      new_positions <- (max_position + 1):(max_position + num_missing_positions)

      # Add all 21 amino acid variants for each new position
      padding_data <- expand.grid(
        mut_aa = all_amino_acids,  # All possible amino acids
        position = new_positions  # New positions to be padded
      )

      # Set placeholder values for the added positions to the exact smallest non-NA value
      padding_data$total_counts_per_cov <- min_non_na_value  # Set to the smallest non-NA value
      padding_data$wt_aa <- "Y"  # Set wild-type amino acid to 'Y'
      padding_data$wt_aa_pos <- paste0("Y", padding_data$position)  # Create wt_aa_pos with correct positions
      padding_data$row_group <- max(heatmap_data_long$row_group)  # Set row group to the current last group

      # Add the new padding rows to heatmap_data_long
      heatmap_data_long <- dplyr::bind_rows(heatmap_data_long, padding_data)
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
  heatmap_data_long <- heatmap_data %>%
    select(mut_aa, position, total_counts_per_cov, wt_aa)  # Use 'total_counts_per_cov'

  # Find the smallest non-NA value in total_counts_per_cov
  min_non_na_value <- min(heatmap_data_long$total_counts_per_cov, na.rm = TRUE)

  # Group positions by rows (75 positions per row) and calculate row_group
  heatmap_data_long <- heatmap_data_long %>%
    mutate(row_group = ((position - 1) %/% 75) + 1)  # Grouping positions into rows

  # Apply padding to add missing positions at the end of the last row, using the calculated min value
  heatmap_data_long <- pad_heatmap_data_long(heatmap_data_long, min_non_na_value)

  # Convert positions to numeric, sort them, and create wt_aa_pos for the plot
  heatmap_data_long <- heatmap_data_long %>%
    mutate(position = as.numeric(position)) %>%  # Ensure position is numeric
    arrange(position) %>%  # Sort by position
    mutate(wt_aa_pos = factor(paste0(wt_aa, position), levels = unique(paste0(wt_aa, position))))  # Create sorted factor levels for wt_aa_pos

  # Add a column to identify synonymous mutations (where mut_aa == wt_aa)
  heatmap_data_long <- heatmap_data_long %>%
    mutate(synonymous = mut_aa == wt_aa)

  # Definiere die korrekte Reihenfolge der Aminosäuren
  amino_acid_order <- c("*", "A", "C", "D", "E", "F", "G", "H", "I", "K",
                        "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

  # Bearbeite heatmap_data_long und erstelle syn_positions gleichzeitig
  syn_positions <- heatmap_data_long %>%
    mutate(mut_aa = factor(mut_aa, levels = amino_acid_order),
           # Berechne die x-Koordinate, die pro Gruppe immer von 1 bis 75 verläuft
           x = as.numeric(factor(wt_aa_pos, levels = unique(wt_aa_pos))) - ((row_group - 1) * 75),
           y = as.numeric(factor(mut_aa, levels = amino_acid_order))) %>%
    filter(synonymous == TRUE)

  # Calculate the number of row groups and adjust plot height dynamically
  num_row_groups <- max(heatmap_data_long$row_group)
  plot_height <- num_row_groups * 4

  # Set the limits for the color scale, ignoring NA (negative values are now NA)
  min_count <- min(heatmap_data_long$total_counts_per_cov, na.rm = TRUE)
  max_count <- max(heatmap_data_long$total_counts_per_cov, na.rm = TRUE)
  max_position <- max(heatmap_data$position)

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
    geom_segment(data = syn_positions[syn_positions$position <= max_position, ],
                 aes(x = x - 0.485, xend = x + 0.485,
                     y = y - 0.485, yend = y + 0.485, color = synonymous),
                 size = 0.2) +

    # Manuelle Farbskala für die diagonalen Linien
    scale_color_manual(values = c("TRUE" = "grey10"), labels = c("TRUE" = "")) +

    theme_minimal() +
    labs(title = "Heatmap of Counts per Coverage for Mutations", x = "Wild-type Amino Acid", y = "Mutant Amino Acid", fill = "Counts per \n Coverage", color = "Synonymous Mutation") +
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

# Test the function
# counts_per_cov_heatmap(
#   "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/prepared_gatk_data.csv",
#   "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/testing_outputs/heatmap.pdf",
#   threshold = 3
# )
