# Load required libraries
library(tidyverse)
library(data.table)
library(qqman)

# Set seed for reproducibility
set.seed(1989)
# ---------------------------
# Define file paths
# ---------------------------
gwas_dir <- "/scratch/ma95362/eth_national_analysis/all_fastq_reads/pangenome_tools_results/bactopia/bactopia-runs/pangenome-20251117-164856/panaroo/filtered_output/pyseer_out"
gwas_file <- file.path(gwas_dir, "New_L4_T3_ETHfamily_gwas.txt")
counts_file <- file.path(gwas_dir, "count_New_L4_T3_ETHfamily_patterns.txt")

cat("Reading GWAS results from:", gwas_file, "\n")
gwas_data <- fread(gwas_file, data.table = FALSE)

# ---------------------------
# Extract number of patterns from counts file
# ---------------------------
counts_lines <- readLines(counts_file)
patterns_line <- counts_lines[grep("^Patterns:", counts_lines)]
num_patterns <- as.numeric(gsub("Patterns:\\s+", "", patterns_line))
cat("Number of patterns detected:", num_patterns, "\n")

# ---------------------------
# Clean and sort results
# ---------------------------
gwas_data <- gwas_data[order(gwas_data$`lrt-pvalue`), ]
gwas_data <- gwas_data[!grepl("bad-chisq", gwas_data$notes), ]

# ---------------------------
# Apply significance threshold (Bonferroni)
# ---------------------------
sig_threshold <- 0.05 / num_patterns
cat("Bonferroni significance threshold:", sig_threshold, "\n")

sig_hits <- gwas_data[gwas_data$`lrt-pvalue` < sig_threshold, ]
cat("Number of significant hits:", nrow(sig_hits), "\n")

# ---------------------------
# Save outputs
# ---------------------------
write.csv(gwas_data,
          file = file.path(gwas_dir, "New_L4_T3_ETHfamily_results.csv"),
          quote = FALSE, row.names = FALSE)

write.csv(sig_hits,
          file = file.path(gwas_dir, "New_L4_T3_ETHfamily_significant_hits.csv"),
          quote = FALSE, row.names = FALSE)

cat("✅ Results saved successfully in", gwas_dir, "\n")