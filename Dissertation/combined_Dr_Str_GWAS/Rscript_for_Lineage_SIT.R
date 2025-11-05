# Load libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(epitools)

# Load metadata
meta <- read_tsv("/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_of_1368/Full_metadata.tab")

# Clean SIT and family columns
meta <- meta %>%
  mutate(
    SIT = gsub("SIT'|'", "", SIT),
    family = gsub("family'|'", "", family)
  )

# Function to run enrichment test and calculate odds ratio
run_enrichment <- function(var, label) {
  df <- meta %>%
    mutate(path1 = Path_ID_cysE_neighbourhood_gene == 1, group = !!sym(var)) %>%
    mutate(group = group == label) %>%
    count(path1, group) %>%
    pivot_wider(names_from = group, values_from = n, values_fill = 0)
  mat <- as.matrix(df[, -1])
  test <- if (any(chisq.test(mat)$expected < 5)) fisher.test(mat) else chisq.test(mat)
  or <- oddsratio(mat)
  data.frame(
    Comparison = paste("Path 1 vs", label),
    Test = ifelse(any(chisq.test(mat)$expected < 5), "Fisher", "Chi-squared"),
    P_value = test$p.value,
    Odds_Ratio = or$measure[2,1],
    CI_lower = or$measure[2,2],
    CI_upper = or$measure[2,3],
    stringsAsFactors = FALSE
  )
}

# Run all comparisons
summary_4222 <- run_enrichment("sub_lineage", "lineage4.2.2.2")
summary_t3eth <- run_enrichment("family", "T3-ETH")
summary_sit149 <- run_enrichment("SIT", "149")

# Combine and save
enrichment_summary <- bind_rows(summary_4222, summary_t3eth, summary_sit149)
write.table(enrichment_summary, "enrichment_summary_with_OR.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Forest plot
ggplot(enrichment_summary, aes(x = Comparison, y = Odds_Ratio)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray") +
  scale_y_log10() +
  labs(title = "Odds Ratios for Path 1 Enrichment", y = "Odds Ratio (log scale)", x = "") +
  theme_minimal() +
  coord_flip()
ggsave("path1_enrichment_forestplot.pdf", width = 8, height = 5)