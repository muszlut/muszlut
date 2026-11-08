#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Annotated phylogeny for the eastern Ethiopia M. bovis / M. caprae dataset.
#
# Input : IQ-TREE tree from the Bactopia snippy run
#         master_sample_table.csv from build_master_table.py
# Output: mbovis_tree.pdf / .png, plus an iTOL annotation file
#
# Run:  Rscript plot_tree.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggtree); library(ggplot2); library(dplyr); library(readr); library(treeio)
})

TREE  <- Sys.glob("/scratch/ma95362/Bio_project_publication/Mbovis/bactopia-runs/snippy-*/core-snp.treefile")[1]
META  <- "/scratch/ma95362/master_sample_table.csv"
OUT   <- "/scratch/ma95362/figures"
dir.create(OUT, showWarnings = FALSE)

stopifnot(file.exists(TREE), file.exists(META))
cat("tree:", TREE, "\n")

tr <- read.tree(TREE)
md <- read_csv(META, show_col_types = FALSE) %>%
  mutate(label = as.character(sample))

# ---- root on the M. caprae clade (the outgroup in this dataset) ------------
cap <- md$label[grepl("caprae|La2", md$rd_species_call, ignore.case = TRUE)]
cap <- intersect(cap, tr$tip.label)
if (length(cap) >= 2) {
  tr <- root(tr, outgroup = cap, resolve.root = TRUE)
  cat("rooted on M. caprae clade:", paste(cap, collapse = ", "), "\n")
} else if ("Reference" %in% tr$tip.label) {
  tr <- root(tr, outgroup = "Reference", resolve.root = TRUE)
}

# ---- bootstrap support -----------------------------------------------------
# IQ-TREE writes UFBoot values as node labels.
supp <- suppressWarnings(as.numeric(tr$node.label))

p <- ggtree(tr, size = 0.4) %<+% md

p <- p +
  # support only where it is meaningful and not maximal
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) &
                              as.numeric(label) >= 95),
                 size = 1.4, colour = "grey25") +
  geom_tippoint(aes(colour = vsnp_clonal_complex), size = 2.2, na.rm = TRUE) +
  geom_tiplab(aes(label = paste0(label, "  ", spoligo_SB)),
              size = 2.3, offset = 0.00002, align = TRUE, linesize = 0.15) +
  scale_colour_manual(
    name   = "Group",
    values = c("Af2"              = "#1f6f8b",
               "non-Af2"          = "#c26a3d",
               "M. caprae (La2)"  = "#6b4c9a",
               "unassigned"       = "grey60"),
    na.value = "grey60") +
  theme_tree2() +
  theme(legend.position = c(0.12, 0.85),
        legend.background = element_rect(fill = "white", colour = "grey80"),
        plot.margin = margin(6, 60, 6, 6)) +
  labs(caption = paste0(
    "Maximum-likelihood phylogeny, IQ-TREE, masked core-SNP alignment. ",
    "Nodes with UFBoot >= 95 marked. Rooted on the M. caprae clade."))

p <- p + ggexpand(0.35, direction = 1)   # room for tip labels

ggsave(file.path(OUT, "mbovis_tree.pdf"), p, width = 9, height = 11)
ggsave(file.path(OUT, "mbovis_tree.png"), p, width = 9, height = 11, dpi = 300)
cat("wrote", file.path(OUT, "mbovis_tree.pdf"), "\n")

# ---- second panel: sublineage + cluster ------------------------------------
clus_col <- grep("^snp_cluster", names(md), value = TRUE)[1]
if (!is.na(clus_col)) {
  hm <- md %>%
    select(label, `Sublineage` = tbp_sublineage,
           `Spoligotype` = spoligo_SB, !!clus_col) %>%
    tibble::column_to_rownames("label")
  p2 <- gheatmap(p, hm, width = 0.25, offset = 0.00012,
                 colnames_angle = 90, colnames_offset_y = 0.4,
                 font.size = 2.5) +
        scale_fill_viridis_d(option = "D", na.value = "white", name = "")
  ggsave(file.path(OUT, "mbovis_tree_annotated.pdf"), p2, width = 12, height = 11)
  cat("wrote", file.path(OUT, "mbovis_tree_annotated.pdf"), "\n")
}

# ---- iTOL colour-strip file ------------------------------------------------
itol <- file.path(OUT, "itol_clonal_complex.txt")
pal <- c("Af2" = "#1f6f8b", "non-Af2" = "#c26a3d",
         "M. caprae (La2)" = "#6b4c9a", "unassigned" = "#999999")
writeLines(c(
  "DATASET_COLORSTRIP", "SEPARATOR TAB", "DATASET_LABEL\tClonal complex",
  "COLOR\t#1f6f8b", "LEGEND_TITLE\tGroup",
  paste0("LEGEND_SHAPES\t", paste(rep(1, length(pal)), collapse = "\t")),
  paste0("LEGEND_COLORS\t", paste(pal, collapse = "\t")),
  paste0("LEGEND_LABELS\t", paste(names(pal), collapse = "\t")),
  "DATA",
  apply(md[!is.na(md$vsnp_clonal_complex), ], 1, function(r)
    paste(r[["label"]], pal[[r[["vsnp_clonal_complex"]]]], sep = "\t"))
), itol)
cat("wrote", itol, "\n")