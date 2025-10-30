# Load libraries
library(tidyverse)
library(ggforce)

# === INPUT ===
summary_counts <- tribble(
  ~Category, ~Count,
  "Core",       3775,
  "Soft core",   84,
  "Shell",      252,
  "Cloud",      591
)

# === CALCULATE PROPORTIONS ===
summary_counts <- summary_counts %>%
  mutate(Percent = Count / sum(Count),
         CumSum = cumsum(Percent))

# === DEFINE LAYER RADII (inner → outer) ===
base_inner <- 0
summary_counts <- summary_counts %>%
  mutate(inner_r = lag(cumsum(Percent), default = base_inner),
         outer_r = cumsum(Percent))

# === ETHIOPIAN FLAG COLORS ===
colors <- c("Core" = "#078930",       # Green
            "Soft core" = "#FFD700",  # Yellow
            "Shell" = "#DA1212",      # Red
            "Cloud" = "#B0B0B0")      # Gray

# === PLOT ===
p <- ggplot(summary_counts) +
  geom_arc_bar(aes(
    x0 = 0, y0 = 0,
    r0 = inner_r, r = outer_r,
    start = 0, end = 2 * pi,
    fill = Category
  ),
  color = "white", size = 0.7
  ) +
  scale_fill_manual(values = colors) +
  coord_fixed() +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  ) +
  labs(
    title = "Ethiopian-Colored Concentric Representation of the Pangenome",
    subtitle = "Innermost = Core → Outermost = Cloud"
  )

# === SAVE OUTPUT ===
ggsave("ethiopian_flag_pangenome.png", plot = p, width = 6, height = 6, dpi = 300)

print(p)
