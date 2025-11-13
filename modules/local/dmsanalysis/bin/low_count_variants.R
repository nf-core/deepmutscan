#!!!UNFINISHED!!! --> Decide if those graphics make sense for users. If so, further develop them and put everything in a function :)


library(tidyverse)

completed_prefiltered_gatk <- read.csv("/Users/benjaminwehnert/CRG/DMS_QC/testing_data/completed_prefiltered_gatk.csv")
View(completed_prefiltered_gatk)

# Lade die notwendigen Pakete
library(dplyr)
library(tidyr)

# Funktion zur Klassifizierung der Basenwechsel
classify_base_changes <- function(wt_codon, variant_codon) {
  # Definition der Purine und Pyrimidine
  purines <- c("A", "G")
  pyrimidines <- c("C", "T")

  # Initialisiere die Zähler
  pur_pur <- 0
  pur_pyr <- 0
  pyr_pur <- 0
  pyr_pyr <- 0

  # Vergleiche jede Base zwischen wt_codon und variant_codon
  for (i in 1:3) {
    wt_base <- substr(wt_codon, i, i)
    var_base <- substr(variant_codon, i, i)

    # Klassifiziere den Basenwechsel
    if (wt_base %in% purines && var_base %in% purines) {
      pur_pur <- pur_pur + 1
    } else if (wt_base %in% purines && var_base %in% pyrimidines) {
      pur_pyr <- pur_pyr + 1
    } else if (wt_base %in% pyrimidines && var_base %in% purines) {
      pyr_pur <- pyr_pur + 1
    } else if (wt_base %in% pyrimidines && var_base %in% pyrimidines) {
      pyr_pyr <- pyr_pyr + 1
    }
  }

  # Erstelle eine Ausgabe als DataFrame
  return(data.frame(pur_pur, pur_pyr, pyr_pur, pyr_pyr))
}

# Beispielanwendung auf den Datensatz
completed_prefiltered_gatk <- completed_prefiltered_gatk %>%
  rowwise() %>%
  mutate(base_change_classification = list(classify_base_changes(wt_codon, Variant))) %>%
  unnest(base_change_classification)








# Load necessary packages
library(dplyr)
library(ggplot2)

# 1. Data preparation
# Create a binary variable for low vs. high counts_per_cov
threshold <- quantile(completed_prefiltered_gatk$counts_per_cov, 0.10, na.rm = TRUE)
completed_prefiltered_gatk <- completed_prefiltered_gatk %>%
  mutate(group = ifelse(counts_per_cov < threshold, "Bottom 10%", "Top 90%"))

# 2. Calculate the mean frequency for each category in both groups
frequency_data <- completed_prefiltered_gatk %>%
  group_by(group) %>%
  summarise(
    pur_pur_mean = mean(pur_pur, na.rm = TRUE),
    pur_pyr_mean = mean(pur_pyr, na.rm = TRUE),
    pyr_pur_mean = mean(pyr_pur, na.rm = TRUE),
    pyr_pyr_mean = mean(pyr_pyr, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = c(pur_pur_mean, pur_pyr_mean, pyr_pur_mean, pyr_pyr_mean),
               names_to = "mutation_type",
               values_to = "mean_frequency") %>%
  # Ensure that all mutation types are present
  complete(mutation_type, group, fill = list(mean_frequency = 0))

# 3. Visualize the distribution of frequencies for the top 90% and bottom 10%
ggplot(frequency_data, aes(x = mutation_type, y = mean_frequency, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Mutation Type", y = "Mean Frequency", title = "Distribution of Mutation Types for Top 90% vs. Bottom 10%") +
  scale_fill_manual(values = c("Top 90%" = "steelblue", "Bottom 10%" = "tomato")) +
  theme_minimal() +
  scale_x_discrete(labels = c("pur_pur" = "Purine > Purine", "pur_pyr" = "Purine > Pyrimidine",
                              "pyr_pur" = "Pyrimidine > Purine", "pyr_pyr" = "Pyrimidine > Pyrimidine"))










threshold <- quantile(completed_prefiltered_gatk$counts_per_cov, 0.10, na.rm = TRUE)
completed_prefiltered_gatk <- completed_prefiltered_gatk %>%
  mutate(group = ifelse(counts_per_cov < threshold, "Bottom 10%", "Top 90%"))

# 2. Create a flag for variants with at least one pur_pyr or pur_pur
completed_prefiltered_gatk <- completed_prefiltered_gatk %>%
  mutate(has_pur_mutation = ifelse(pyr_pur > 0 | pur_pyr > 0, 1, 0))

# 3. Calculate the proportion of variants with at least one pur_pyr or pur_pur for each group
proportion_data <- completed_prefiltered_gatk %>%
  group_by(group) %>%
  summarise(
    total_variants = n(),
    variants_with_pur_mutation = sum(has_pur_mutation),
    proportion_with_pur_mutation = variants_with_pur_mutation / total_variants
  )

# 4. Print the results
print(proportion_data)

# 5. Optional: Visualize the proportion for each group
ggplot(proportion_data, aes(x = group, y = proportion_with_pur_mutation, fill = group)) +
  geom_bar(stat = "identity") +
  labs(x = "Group", y = "Proportion of Variants", title = "Proportion of Variants with at Least One Pur_pyr or Pur_pur Mutation") +
  scale_fill_manual(values = c("Top 90%" = "steelblue", "Bottom 10%" = "tomato")) +
  theme_minimal()



