# Load libraries
library(dplyr)
library(ggplot2)
library(readr)

# Load metadata
meta <- read_tsv("/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_of_1368/Full_metadata.tab")


# Count genomes per path and sublineage
path_lineage_counts <- meta %>%
  count(Path_ID_cysE_neighbourhood_gene, sub_lineage) %>%
  rename(Path_ID = Path_ID_cysE_neighbourhood_gene)

# Plot stacked barplot: path vs sublineage
ggplot(path_lineage_counts, aes(x = factor(Path_ID), y = n, fill = sub_lineage)) +
  geom_bar(stat = "identity") +
  labs(title = "Path–Sublineage Association",
       x = "Path ID (cysE neighborhood)",
       y = "Number of Genomes",
       fill = "Sublineage") +
  theme_minimal()

# Optional: test enrichment of Path 1 in sublineage 4.2.2.2
enrichment_table <- meta %>%
  mutate(path1 = Path_ID_cysE_neighbourhood_gene == 1,
         is_4222 = sub_lineage == "lineage4.2.2.2") %>%
  count(path1, is_4222)

chisq.test(enrichment_table)