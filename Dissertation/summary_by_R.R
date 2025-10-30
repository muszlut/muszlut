# Load necessary package
library(tidyverse)

# === INPUT FILES ===
# Gene presence/absence file from Panaroo
# Typically named something like: "gene_presence_absence.Rtab"
rtab_file <- "gene_presence_absence.Rtab"

# === PARAMETERS ===
# Number of strains in your dataset (can be counted from Rtab)
# This is important to compute gene frequency correctly.
num_strains <- length(read_tsv(rtab_file, n_max = 0)) - 1  # subtract the 'Gene' column

# === LOAD DATA ===
data <- read_tsv(rtab_file)

# === CALCULATE GENE PRESENCE FREQUENCY ===
# For each gene, calculate the number of strains where the gene is present (value == 1)
data$Presence_Count <- rowSums(data[,-1] != 0)  # assuming 1=present, 0=absent
data$Presence_Percent <- (data$Presence_Count / num_strains) * 100

# === CATEGORIZE GENES BASED ON YOUR SUMMARY CRITERIA ===
data$Category <- case_when(
  data$Presence_Percent >= 99 ~ "Core genes (99-100%)",
  data$Presence_Percent >= 95 & data$Presence_Percent < 99 ~ "Soft core (95-99%)",
  data$Presence_Percent >= 15 & data$Presence_Percent < 95 ~ "Shell (15-95%)",
  data$Presence_Percent < 15 ~ "Cloud (<15%)"
)

# === SAVE RESULTS ===
write_tsv(data, "gene_categorization_panaroo.tsv")

# === OPTIONAL: SUMMARY COUNTS ===
summary_counts <- data %>%
  group_by(Category) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))

print(summary_counts)
write_tsv(summary_counts, "gene_category_summary.tsv")

# === OPTIONAL: VISUALIZE DISTRIBUTION ===
library(ggplot2)
ggplot(summary_counts, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  coord_flip() +
  labs(title = "Panaroo Gene Categorization Summary",
       x = "Category",
       y = "Number of Genes")
