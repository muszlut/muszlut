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

# Plot: Path vs Sublineage
p1 <- ggplot(path_lineage_counts, aes(x = factor(Path_ID), y = n, fill = sub_lineage)) +
  geom_bar(stat = "identity") +
  labs(title = "Path - Sublineage Association", x = "Path ID", y = "Number of Genomes", fill = "Sublineage") +
  theme_minimal()
ggsave("path_sublineage_plot.pdf", plot = p1, width = 8, height = 6)

# Highlight SIT-149 or T3-ETH
highlight_group <- meta %>% filter(SIT == "149" | family == "T3-ETH")

# Plot by Path ID and Phenotype
p2 <- ggplot(highlight_group, aes(x = factor(Path_ID_cysE_neighbourhood_gene), fill = Phenotype)) +
  geom_bar(position = "stack") +
  labs(title = "SIT-149 / T3-ETH by Path ID and Phenotype", x = "Path ID", y = "Number of Genomes", fill = "Phenotype") +
  theme_minimal()
ggsave("SIT149_T3ETH_by_Phenotype.pdf", plot = p2, width = 8, height = 6)

# Plot by Path ID and Drug Resistance Class
p3 <- ggplot(highlight_group, aes(x = factor(Path_ID_cysE_neighbourhood_gene), fill = dr_class)) +
  geom_bar(position = "stack") +
  labs(title = "SIT-149 / T3-ETH by Path ID and Drug Resistance", x = "Path ID", y = "Number of Genomes", fill = "Drug Resistance Class") +
  theme_minimal()
ggsave("SIT149_T3ETH_by_dr_class.pdf", plot = p3, width = 8, height = 6)

# Enrichment test: Path 1 vs sublineage 4.2.2.2
enrichment_table <- meta %>%
  mutate(
    path1 = Path_ID_cysE_neighbourhood_gene == 1,
    is_4222 = sub_lineage == "lineage4.2.2.2"
  ) %>%
  count(path1, is_4222) %>%
  pivot_wider(names_from = is_4222, values_from = n, values_fill = 0)

test_matrix <- as.matrix(enrichment_table[, -1])
chi <- chisq.test(test_matrix)

if (any(chi$expected < 5)) {
  message("Using Fisher's exact test due to low expected counts")
  fisher_result <- fisher.test(test_matrix)
  capture.output(fisher_result, file = "fisher_result.txt")
} else {
  message("Using Chi-squared test")
  capture.output(chi, file = "chi_squared_result.txt")
}

# Save enrichment table
write.table(test_matrix, "enrichment_table.txt", sep = "\t", quote = FALSE)