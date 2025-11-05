# Load required libraries
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)

# Load metadata
meta <- read_tsv("/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_of_1368/Full_metadata.tab")

# Clean SIT and family columns
meta <- meta %>%
  mutate(
    SIT = gsub("SIT'|'", "", SIT),
    family = gsub("family'|'", "", family)
  )

# Count genomes per path and sublineage
path_lineage_counts <- meta %>%
  count(Path_ID_cysE_neighbourhood_gene, sub_lineage) %>%
  rename(Path_ID = Path_ID_cysE_neighbourhood_gene)

# Stacked barplot: Path ID vs Sublineage
ggplot(path_lineage_counts, aes(x = factor(Path_ID), y = n, fill = sub_lineage)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Path - Sublineage Association",
    x = "Path ID (cysE neighborhood)",
    y = "Number of Genomes",
    fill = "Sublineage"
  ) +
  theme_minimal()

# Filter for SIT-149 or family T3-ETH
highlight_group <- meta %>%
  filter(SIT == "149" | family == "T3-ETH")

# Barplot: Highlighted group by Path ID and Phenotype
ggplot(highlight_group, aes(x = factor(Path_ID_cysE_neighbourhood_gene), fill = Phenotype)) +
  geom_bar(position = "stack") +
  labs(
    title = "SIT-149 / T3-ETH Group by Path ID",
    x = "Path ID",
    y = "Number of Genomes",
    fill = "Phenotype"
  ) +
  theme_minimal()

# Enrichment test: Path 1 vs sublineage 4.2.2.2
enrichment_table <- meta %>%
  mutate(
    path1 = Path_ID_cysE_neighbourhood_gene == 1,
    is_4222 = sub_lineage == "lineage4.2.2.2"
  ) %>%
  count(path1, is_4222) %>%
  pivot_wider(names_from = is_4222, values_from = n, values_fill = 0)

# Convert to matrix and run Chi-squared test
test_matrix <- as.matrix(enrichment_table[, -1])
chi <- chisq.test(test_matrix)

# Use Fisher's exact test if expected counts are low
if (any(chi$expected < 5)) {
  message("Using Fisher's exact test due to low expected counts")
  fisher.test(test_matrix)
} else {
  message("Using Chi-squared test")
  chi
}