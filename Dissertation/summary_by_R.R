# Load libraries
library(tidyverse)
library(ggplot2)

# === INPUT FILE ===
rtab_file <- "gene_presence_absence.Rtab"

# === READ AND PROCESS DATA ===
num_strains <- length(read_tsv(rtab_file, n_max = 0)) - 1
data <- read_tsv(rtab_file)

# Calculate presence %
data$Presence_Count <- rowSums(data[,-1] != 0)
data$Presence_Percent <- (data$Presence_Count / num_strains) * 100

# Categorize based on Panaroo thresholds
data$Category <- case_when(
  data$Presence_Percent >= 99 ~ "Core",
  data$Presence_Percent >= 95 & data$Presence_Percent < 99 ~ "Soft core",
  data$Presence_Percent >= 15 & data$Presence_Percent < 95 ~ "Shell",
  data$Presence_Percent < 15 ~ "Cloud"
)

# === SUMMARIZE CATEGORIES ===
summary_counts <- data %>%
  group_by(Category) %>%
  summarise(Count = n()) %>%
  mutate(Percent = Count / sum(Count) * 100) %>%
  arrange(match(Category, c("Cloud", "Shell", "Soft core", "Core")))

print(summary_counts)

# === CREATE PROPORTIONAL CIRCULAR LAYERS ===
# Calculate ring radii based on cumulative proportions
summary_counts <- summary_counts %>%
  mutate(
    outer_radius = cumsum(Percent),
    inner_radius = lag(outer_radius, default = 0)
  )

# Plot using ggplot2 with geom_rect + coord_polar
ggplot() +
  geom_rect(data = summary_counts,
            aes(xmin = 0, xmax = 1,
                ymin = inner_radius, ymax = outer_radius,
                fill = Category),
            color = "white", linewidth = 0.5) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c(
    "Cloud" = "#ffb3b3",      # outermost
    "Shell" = "#ffcc66",      # orange-yellow
    "Soft core" = "#99cc99",  # green
    "Core" = "#3399ff"        # innermost
  )) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  ) +
  labs(title = "Proportional Circular Pangenome Composition",
       subtitle = "Innermost = Core → Outermost = Cloud") +
  annotate("text", x = 0, y = mean(summary_counts$outer_radius), label = "")

# === SAVE OUTPUTS ===
write_tsv(data, "gene_categorization_panaroo.tsv")
write_tsv(summary_counts, "gene_category_summary.tsv")
ggsave("proportional_circular_pangenome_layers.png", width = 6, height = 6, dpi = 300)
cat("✅ Pangenome summary and plot saved successfully.\n")