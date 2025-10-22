library(tidyverse)


####set seed to make this process more reproducible####


set.seed(1989)

library(data.table)
library(qqman)

# Read GWAS result
tbl <- fread("T3-ETH_gwas.txt", data.table = FALSE)

# Remove bad-chisq entries
tbl <- tbl[!grepl("bad-chisq", tbl$notes), ]

# Set threshold
sig_threshold <- 1.21e-4

# Filter significant hits
sig_hits <- subset(tbl, `lrt-pvalue` < sig_threshold)

# Save results
write.csv(sig_hits, "T3-ETH_significant_hitsbyR.csv", row.names = FALSE, quote = FALSE)

# Optional: Manhattan and QQ plots
manhattan(tbl, p="lrt-pvalue", main="T3-ETH vs others (Pyseer GWAS)")
qq(tbl$`lrt-pvalue`)
qq(sig_hits$`lrt-pvalue`)   