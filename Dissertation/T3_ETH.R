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

# Add dummy CHR and BP columns for plotting
tbl$CHR <- 1  # Assign all to chromosome 1
tbl$BP <- seq_along(tbl$variant)  # Use row index as base pair position
tbl$SNP <- tbl$variant  # Optional: label with variant name

# Manhattan plot with simulated coordinates
manhattan(tbl, chr = "CHR", bp = "BP", snp = "SNP", p = "lrt-pvalue",
          main = "T3ETH vs others (Pyseer GWAS)", col = c("blue4", "orange3"))

# Generate QQ plots
qq(tbl$`lrt-pvalue`, main = "QQ Plot: All Variants")
qq(sig_hits$`lrt-pvalue`, main = "QQ Plot: Significant Variants")