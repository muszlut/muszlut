# ✅ Pre-flight package check
required <- c("ggplot2", "ggrepel", "dplyr", "tidyr", "readr")
missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse = ", "))
}

# 📦 Load libraries
library(ggplot2)
library(dplyr)
library(readr)
library(ggrepel)

# 🔧 Increase buffer size
Sys.setenv("VROOM_CONNECTION_SIZE" = 5000000)

# 🧪 Validate input files
stopifnot(file.exists("alignment_entropy.csv"))
stopifnot(file.exists("core-genome.position_cross_reference.txt.gz"))
stopifnot(file.exists("core-genome.importation_status.txt"))
stopifnot(file.exists("L4.2.2.2_Binary_T3_ETHfamily_significant_hits.csv"))

# ✅ Now load entropy safely
entropy <- read_csv("alignment_entropy.csv", col_names = c("gene", "entropy"), show_col_types = FALSE)

# 2️⃣ Load gene coordinates
coords <- read_tsv("core-genome.position_cross_reference.txt.gz", show_col_types = FALSE)

# 🔗 Merge entropy with coordinates
genes <- coords %>%
  left_join(entropy, by = "gene") %>%
  mutate(entropy_color = case_when(
    is.na(entropy) ~ "Missing",
    entropy > quantile(entropy, 0.9, na.rm = TRUE) ~ "High",
    TRUE ~ "Low"
  ))

# 3️⃣ Load recombinant tracts
recomb <- read_tsv("core-genome.importation_status.txt", show_col_types = FALSE)
recomb$color <- "Recombination"

# 4️⃣ Load pyseer hits
pyseer <- read_tsv("L4.2.2.2_Binary_T3_ETHfamily_significant_hits.csv", show_col_types = FALSE) %>%
  select(gene, beta, lrt_pvalue)

# 🔗 Merge pyseer hits with gene coordinates
pyseer <- coords %>%
  inner_join(pyseer, by = "gene") %>%
  mutate(y_pos = 1.05)

# 5️⃣ Base plot
p <- ggplot() +
  geom_segment(data = genes, aes(x = start, xend = end, y = 0, yend = 0, color = entropy_color), size = 4) +
  scale_color_manual(values = c("High" = "red", "Low" = "grey", "Missing" = "lightblue")) +
  geom_segment(data = recomb, aes(x = Beg, xend = End, y = 0.1, yend = 0.1), color = "blue", size = 2, alpha = 0.7) +
  geom_point(data = pyseer, aes(x = start, y = y_pos), shape = 8, color = "black", size = 3) +
  geom_text_repel(data = genes %>% filter(entropy > quantile(entropy, 0.95, na.rm = TRUE)),
                  aes(x = (start + end) / 2, y = 0.05, label = gene), size = 3) +
  theme_bw() +
  labs(x = "Genome position (bp)", y = "", color = "Entropy") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# 6️⃣ Save plot
ggsave("Genome_entropy_recomb_pyseer.png", p, width = 12, height = 3, dpi = 300)

# 📋 Save session info
writeLines(capture.output(sessionInfo()), "session_info.txt")