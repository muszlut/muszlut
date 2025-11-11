# ===============================
# R Script: Genome-wide plot
# ===============================
library(ggplot2)
library(dplyr)
library(readr)
library(ggrepel)

# 1️⃣ Load gene-level entropy
entropy <- read_tsv("alignment_entropy.csv") 

# 2️⃣ Load gene coordinates
coords <- read_tsv("core-genome.position_cross_reference.txt.gz") 

# Merge entropy with coordinates
genes <- coords %>%
  left_join(entropy, by="gene") %>%
  mutate(entropy_color = ifelse(entropy > quantile(entropy,0.9, na.rm=TRUE), "High", "Low"))

# 3️⃣ Load recombinant tracts
recomb <- read_tsv("core-genome.importation_status.txt") 

# Optional: assign a single color for recombination
recomb$color <- "Recombination"

# 4️⃣ Load pyseer hits
pyseer <- read_tsv("L4.2.2.2_Binary_T3_ETHfamily_significant_hits.csv") %>%
  select(gene, beta, lrt_pvalue)

# Merge pyseer hits with gene coordinates
pyseer <- coords %>% 
  inner_join(pyseer, by="gene") %>%
  mutate(y_pos = 1.05)  # vertical position for star

# 5️⃣ Base plot
p <- ggplot() + 
  # Genes as segments
  geom_segment(data=genes, aes(x=start, xend=end, y=0, yend=0, color=entropy_color), size=4) +
  scale_color_manual(values=c("High"="red","Low"="grey")) +
  
  # Recombinant tracts as points
  geom_segment(data=recomb, aes(x=Beg, xend=End, y=0.1, yend=0.1), color="blue", size=2, alpha=0.7) +
  
  # Pyseer hits as stars
  geom_point(data=pyseer, aes(x=start, y=y_pos), shape=8, color="black", size=3) +
  
  # Labels for top high-entropy genes
  geom_text_repel(data=genes %>% filter(entropy > quantile(entropy,0.95, na.rm=TRUE)),
                  aes(x=(start+end)/2, y=0.05, label=gene), size=3) +
  
  theme_bw() +
  labs(x="Genome position (bp)", y="", color="Entropy") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# 6️⃣ Save plot
ggsave("Genome_entropy_recomb_pyseer.png", p, width=12, height=3, dpi=300)

# ===============================
