# Load required libraries
library(tidyverse)
library(data.table)
library(qqman)

# Set seed for reproducibility
set.seed(1989)

# Read GWAS results
tbl <- fread("T3ETH_gwas.txt", data.table = FALSE)

# Remove entries flagged as bad-chisq
tbl <- tbl[!grepl("bad-chisq", tbl$notes), ]

# Define Bonferroni-corrected significance threshold
sig_threshold <- 1.21e-4

# Filter for significant hits
sig_hits <- subset(tbl, `lrt-pvalue` < sig_threshold)

# Save filtered results
write.csv(sig_hits, "T3ETH_significant_hitsbyR.csv", row.names = FALSE, quote = FALSE)

# Generate Manhattan plot
manhattan(tbl, p = "lrt-pvalue", main = "T3ETH vs others (Pyseer GWAS)")

# Generate QQ plots
qq(tbl$`lrt-pvalue`, main = "QQ Plot: All Variants")
qq(sig_hits$`lrt-pvalue`, main = "QQ Plot: Significant Variants")