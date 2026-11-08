#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Annotated phylogeny using only ape + base graphics.
#
# Deliberately dependency-light: ape is the single requirement, so this runs
# anywhere R does. Produces the same information as a ggtree figure (clonal
# complex by colour, spoligotype on the tip label, bootstrap support on
# well-supported nodes) without the ggtree dependency chain, which cannot be
# satisfied in this R installation (ggtree requires a check_linewidth()
# that no released ggfun provides).
#
# Input : IQ-TREE tree from the Bactopia snippy run. No metadata file is
#         needed -- group membership and spoligotypes are recorded below,
#         taken from the vSNP3 step 2 defining-SNP groups and step 1 stats.
# Output: mbovis_tree_ape.pdf / .png, and an iTOL colour-strip file.
#
# Run:  Rscript --vanilla plot_tree_ape.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(ape))

TREE <- Sys.glob("/scratch/ma95362/Bio_project_publication/Mbovis/bactopia-runs/snippy-*/core-snp.treefile")[1]
OUT  <- "/scratch/ma95362/Bio_project_publication/Mbovis/R_output/figures"

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(TREE))
cat("tree  :", TREE, "\n")
cat("output:", OUT, "\n\n")

# ---- group membership (vSNP3 step 2 defining-SNP groups) ------------------
af2 <- c("17","151P","213P","262G","264G","264P","271G","281G","281P","283G",
         "292G","300P","301G","304G","304P","305G","308P","310P","342P","343P")
non <- c("120PG","231G","258G","259G","262P","288P","307G","318G","322G",
         "323G","32_HK","341P")
cap <- c("097G_400","100GP","106GP","97G")
una <- c("147","234GP")

# ---- spoligotypes (vSNP3 step 1 stats) ------------------------------------
sb <- c(
  SB0134 = "120PG,231G,258G,259G,262P,288P,307G,318G,322G,323G,32_HK,341P",
  SB1176 = "151P,17,262G,264G,264P,271G,292G,300P,301G,304G,304P,308P,310P,342P,343P",
  SB0133 = "281G,281P,283G",
  SB0418 = "097G_400,100GP,106GP",
  SB1690 = "97G"
)
spoligo <- setNames(rep(names(sb), lengths(strsplit(sb, ","))),
                    unlist(strsplit(sb, ",")))

tr <- read.tree(TREE)
cat("tips in tree:", length(tr$tip.label), "\n")

# ---- root on the M. caprae clade ------------------------------------------
og <- intersect(cap, tr$tip.label)
if (length(og) >= 2) {
  tr <- root(tr, outgroup = og, resolve.root = TRUE)
  cat("rooted on M. caprae:", paste(og, collapse = ", "), "\n")
} else if ("Reference" %in% tr$tip.label) {
  tr <- root(tr, outgroup = "Reference", resolve.root = TRUE)
  cat("rooted on Reference\n")
}
tr <- ladderize(tr)

# ---- colours ---------------------------------------------------------------
pal <- c(Af2         = "#1f6f8b",
         `non-Af2`   = "#c26a3d",
         `M. caprae` = "#6b4c9a",
         unassigned  = "#999999",
         Reference   = "#000000")

grp <- setNames(rep("unassigned", length(tr$tip.label)), tr$tip.label)
grp[tr$tip.label %in% af2]      <- "Af2"
grp[tr$tip.label %in% non]      <- "non-Af2"
grp[tr$tip.label %in% cap]      <- "M. caprae"
grp[tr$tip.label == "Reference"] <- "Reference"
tipcol <- pal[grp]

cat("\ngroup counts:\n"); print(table(grp))

# any tip not accounted for?
unknown <- setdiff(tr$tip.label, c(af2, non, cap, una, "Reference"))
if (length(unknown)) cat("\nNOT IN ANY GROUP:", paste(unknown, collapse = ", "), "\n")

# ---- tip labels with spoligotype -------------------------------------------
orig <- tr$tip.label
sbv  <- spoligo[orig]
tr$tip.label <- ifelse(is.na(sbv), orig, sprintf("%-9s %s", orig, sbv))

# ---- draw ------------------------------------------------------------------
for (dev in c("pdf", "png")) {
  f <- file.path(OUT, paste0("mbovis_tree_ape.", dev))
  if (dev == "pdf") pdf(f, width = 9, height = 11)
  else              png(f, width = 9, height = 11, units = "in", res = 300)

  par(mar = c(5, 1, 3, 1), xpd = NA)
  plot.phylo(tr, cex = 0.72, font = 1, label.offset = 2,
             tip.color = tipcol, edge.width = 1.1)

  bs <- suppressWarnings(as.numeric(tr$node.label))
  ok <- !is.na(bs) & bs >= 95
  if (any(ok)) {
    nodes <- (Ntip(tr) + 1):(Ntip(tr) + tr$Nnode)
    nodelabels(node = nodes[ok], pch = 16, cex = 0.5, col = "grey25")
  }

  add.scale.bar(lwd = 1.2, cex = 0.8)
  legend("bottomleft", bty = "n", cex = 0.8, pch = 15, col = pal,
         legend = names(pal), title = "Group")
  legend("bottomright", bty = "n", cex = 0.7, pch = 16, col = "grey25",
         legend = expression("UFBoot" >= 95))
  title(main = expression(italic("M. bovis")*" / "*italic("M. caprae")*", eastern Ethiopia"),
        cex.main = 1.0, font.main = 1)
  mtext(paste("ML phylogeny (IQ-TREE 3.1.3) from the core-SNP alignment;",
              "rooted on the M. caprae clade"),
        side = 1, line = 3.2, cex = 0.7)
  invisible(dev.off())
  cat("wrote", f, "\n")
}

# ---- iTOL colour strip -----------------------------------------------------
itol <- file.path(OUT, "itol_clonal_complex.txt")
ipal <- pal[c("Af2", "non-Af2", "M. caprae", "unassigned")]
writeLines(c(
  "DATASET_COLORSTRIP", "SEPARATOR TAB",
  "DATASET_LABEL\tClonal complex", "COLOR\t#1f6f8b",
  "LEGEND_TITLE\tGroup",
  paste0("LEGEND_SHAPES\t", paste(rep(1, length(ipal)), collapse = "\t")),
  paste0("LEGEND_COLORS\t", paste(ipal, collapse = "\t")),
  paste0("LEGEND_LABELS\t", paste(names(ipal), collapse = "\t")),
  "DATA",
  vapply(orig[orig != "Reference"],
         function(s) paste(s, pal[[grp[[s]]]], sep = "\t"), character(1))
), itol)
cat("wrote", itol, "\n")