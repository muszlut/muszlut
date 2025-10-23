# Load required libraries
library(tidyverse)
library(data.table)
library(qqman)

# Set seed for reproducibility
set.seed(1989)

# Read GWAS results
tbl <- fread("resistance_level_gwas.txt", data.table = FALSE)
tbl <- tbl[!grepl("bad-chisq", tbl$notes), ]

# Define Bonferroni threshold
sig_threshold <- 1.21e-4
sig_hits <- subset(tbl, `lrt-pvalue` < sig_threshold)

# Save filtered results
write.csv(sig_hits, "resistance_level_gwas_significant_hitsbyR.csv", row.names = FALSE, quote = FALSE)

# Add dummy CHR and BP columns for plotting
tbl$CHR <- 1
tbl$BP <- seq_along(tbl$variant)
tbl$SNP <- tbl$variant

# Highlight top hits (e.g., mmaA4 if present)
#highlight_variants <- c("mmaA4")
#highlight_variants <- intersect(highlight_variants, tbl$SNP)

# Save Manhattan plot
png("DR_manhattan_plot.png", width = 1200, height = 800)
manhattan(tbl,
          chr = "CHR", bp = "BP", snp = "SNP", p = "lrt-pvalue",
          main = "DR vs Others (Pyseer GWAS)",
          col = c("steelblue3", "darkorange2"),
          cex = 1.2, cex.axis = 1.2, cex.lab = 1.4,
          suggestiveline = FALSE,
          genomewideline = -log10(sig_threshold),
          #highlight = highlight_variants)
dev.off()

# Save QQ plot (all variants)
png("DR_qq_all.png", width = 800, height = 800)
qq(tbl$`lrt-pvalue`, main = "QQ Plot: All Variants", cex = 1.2)
dev.off()

# Save QQ plot (significant hits)
png("DR_qq_significant.png", width = 800, height = 800)
qq(sig_hits$`lrt-pvalue`, main = "QQ Plot: Significant Variants", cex = 1.2)
dev.off()