# ============================================================================
# SNP-distance clustering, cluster monophyly testing, and pairwise SNP
# distance distribution.
#
# Abdela et al., Microbial Genomics (MGEN-D-26-00378), revision 1.
#
# Inputs (from the Bactopia snippy workflow output directory):
#   core-snp.masked.distance.tsv   snp-dists matrix from the masked alignment
#   core-snp.treefile              IQ-TREE maximum-likelihood tree
#   sample_map.csv                 two columns: srr, sample   (see note below)
#   Book1.xlsx                     isolate metadata (optional, for comparison)
#
# Outputs:
#   cluster_assignments.csv        per-isolate cluster membership
#   cluster_summary.csv            per-cluster size, lineage, monophyly
#   clustering_thresholds.csv      clustered isolates / clusters / RTI by cut-off
#   snp_distance_histogram.png     pairwise SNP distance distribution
#
# NOTE ON SAMPLE NAMING
# The distance matrix and tree use SRA run accessions (SRR...), while the
# metadata uses internal sample IDs (P111, E24, ...). A mapping file is
# required. If you do not have one, it can be built from the Bactopia
# samples file or the SRA submission metadata.
# ============================================================================

library(tidyverse)
library(igraph)
library(ape)
library(readxl)

THRESHOLDS <- c(4, 5, 10, 12)   # 10 = primary; 4 = sensitivity analysis
PRIMARY    <- 10

# ---------------------------------------------------------------------------
# 1. Read the SNP distance matrix
# ---------------------------------------------------------------------------
# snp-dists writes a square matrix whose first column holds the row names and
# whose first header cell contains the program name rather than a column name.

read_snpdists <- function(path) {
  m <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
                  row.names = 1, comment.char = "")
  m <- as.matrix(m)
  storage.mode(m) <- "numeric"
  stopifnot(nrow(m) == ncol(m))
  m
}

dm <- read_snpdists("core-snp.masked.distance.tsv")
message("Matrix read: ", nrow(dm), " sequences")

# Drop the H37Rv reference, which snippy-core includes as a row/column.
# It must not participate in clustering between isolates.
ref_rows <- grep("^Reference$|GCF_|H37Rv", rownames(dm), ignore.case = TRUE)
if (length(ref_rows) > 0) {
  message("Removing reference row(s): ", paste(rownames(dm)[ref_rows], collapse = ", "))
  dm <- dm[-ref_rows, -ref_rows, drop = FALSE]
}
message("Isolates retained for clustering: ", nrow(dm))   # expect 257

# ---------------------------------------------------------------------------
# 2. Single-linkage clustering
# ---------------------------------------------------------------------------
# Isolates separated by <= threshold SNPs are linked; clusters are the
# connected components of the resulting graph. Isolates with no neighbour
# within the threshold are unclustered (NOT singleton clusters).

cluster_at <- function(dm, threshold) {
  adj <- dm <= threshold
  diag(adj) <- FALSE
  g    <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  comp <- components(g)

  memb  <- comp$membership
  sizes <- comp$csize[memb]

  out <- tibble(id = rownames(dm), raw = memb, size = sizes) |>
    mutate(cluster = if_else(size >= 2, raw, NA_integer_))

  # Renumber clusters by size (largest = 1), ties broken by smallest member,
  # so that cluster IDs are stable and reproducible across reruns.
  key <- out |>
    filter(!is.na(cluster)) |>
    group_by(cluster) |>
    summarise(n = n(), first_id = min(id), .groups = "drop") |>
    arrange(desc(n), first_id) |>
    mutate(new = row_number())

  out |>
    left_join(key, by = "cluster") |>
    transmute(id, cluster = new, cluster_size = n)
}

summarise_clustering <- function(cl, n_total) {
  n_clustered <- sum(!is.na(cl$cluster))
  n_clusters  <- n_distinct(cl$cluster[!is.na(cl$cluster)])
  tibble(
    clustered_isolates = n_clustered,
    pct_clustered      = round(100 * n_clustered / n_total, 1),
    n_clusters         = n_clusters,
    RTI                = round(100 * (n_clustered - n_clusters) / n_total, 1)
  )
}

thresh_tbl <- map_dfr(THRESHOLDS, function(th) {
  bind_cols(tibble(threshold_snps = th),
            summarise_clustering(cluster_at(dm, th), nrow(dm)))
})
print(as.data.frame(thresh_tbl), row.names = FALSE)
write_csv(thresh_tbl, "clustering_thresholds.csv")

# EXPECTED (verify against the published Table 1):
#   threshold 10 -> 125 clustered (48.6%), 39 clusters, RTI 33.5%
#   threshold  4 -> 111 clustered (43.2%), 39 clusters, RTI 28.0%

clusters <- cluster_at(dm, PRIMARY)

# ---------------------------------------------------------------------------
# 3. Attach metadata and compare with the previously published assignment
# ---------------------------------------------------------------------------

if (file.exists("Supplementary_Table_S1_accessions.xlsx") && file.exists("Book1.xlsx")) {

  map <- read_excel("Supplementary_Table_S1_accessions.xlsx") |>
  transmute(srr = SRA_accession,
            sample = str_remove(str_trim(library_ID), "\\.$"))   # srr, sample
  meta <- read_excel("Book1.xlsx") |>
    mutate(sample = str_remove(sample, "\\.$"))                # strip trailing dots

  clusters <- clusters |>
    left_join(map, by = c("id" = "srr")) |>
    left_join(meta, by = "sample")

  # Concordance with the web-tool clustering recorded in Book1.xlsx.
  # Cluster NUMBERS will differ (numbering here is by size); what must agree
  # is the PARTITION, i.e. which isolates are grouped together.
  old <- clusters |>
    mutate(old = if_else(str_detect(tolower(Cluster_10), "na"), NA_character_, Cluster_10))

  agree <- old |>
    filter(!is.na(cluster), !is.na(old)) |>
    group_by(cluster) |>
    summarise(n_old_groups = n_distinct(old), .groups = "drop")

  message("\nClusters here that merge >1 previous cluster: ",
          sum(agree$n_old_groups > 1))
  message("Isolates clustered here but not previously: ",
          sum(!is.na(old$cluster) & is.na(old$old)))
  message("Isolates clustered previously but not here: ",
          sum(is.na(old$cluster) & !is.na(old$old)))
}

write_csv(clusters, "cluster_assignments.csv")

# ---------------------------------------------------------------------------
# 4. Monophyly of each cluster on the ML tree
# ---------------------------------------------------------------------------
# This identifies the polyphyletic cluster the Editor and both Reviewers asked
# about. A cluster is monophyletic if its members form a clade containing no
# other isolate. Single-linkage clustering does not guarantee this: chains of
# pairwise links can group isolates that are not reciprocally within the
# threshold.

tr <- read.tree("core-snp.treefile")
tr <- drop.tip(tr, grep("^Reference$|GCF_|H37Rv", tr$tip.label,
                        ignore.case = TRUE, value = TRUE))

mono <- clusters |>
  filter(!is.na(cluster)) |>
  group_by(cluster) |>
  summarise(n = n(), tips = list(id), .groups = "drop") |>
  mutate(
    monophyletic = map_lgl(tips, ~ {
      tp <- intersect(.x, tr$tip.label)
      if (length(tp) < 2 || length(tp) < length(.x)) return(NA)
      is.monophyletic(tr, tp)
    }),
    members = map_chr(tips, ~ paste(sort(.x), collapse = "; "))
  ) |>
  select(-tips)

message("\n--- Cluster monophyly ---")
print(as.data.frame(select(mono, cluster, n, monophyletic)), row.names = FALSE)
message("\nPOLYPHYLETIC CLUSTER(S): ",
        paste(mono$cluster[which(!mono$monophyletic)], collapse = ", "))

# For each polyphyletic cluster, list the intruding isolates: tips that fall
# inside the clade spanning the cluster but are not cluster members.
for (cc in mono$cluster[which(!mono$monophyletic)]) {
  tp   <- clusters$id[which(clusters$cluster == cc)]
  node <- getMRCA(tr, intersect(tp, tr$tip.label))
  clade <- extract.clade(tr, node)$tip.label
  message("Cluster ", cc, " MRCA clade contains ", length(clade),
          " tips; non-members: ", paste(setdiff(clade, tp), collapse = ", "))
}

write_csv(mono, "cluster_summary.csv")

# ---------------------------------------------------------------------------
# 5. Pairwise SNP distance distribution  (requested by Reviewer 2)
# ---------------------------------------------------------------------------

d <- dm[upper.tri(dm)]
message("\nPairwise distances: n = ", length(d),
        "; median = ", median(d), "; range ", min(d), "-", max(d))
message("Pairs <= 10 SNPs: ", sum(d <= 10),
        " (", round(100 * sum(d <= 10) / length(d), 3), "%)")

p_all <- ggplot(tibble(d = d), aes(d)) +
  geom_histogram(binwidth = 25, fill = "grey40", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = PRIMARY, colour = "firebrick", linetype = "dashed") +
  labs(x = "Pairwise SNP distance", y = "Number of isolate pairs",
       title = "All pairwise distances") +
  theme_classic(base_size = 11)

p_zoom <- ggplot(tibble(d = d[d <= 100]), aes(d)) +
  geom_histogram(binwidth = 1, fill = "steelblue", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = PRIMARY, colour = "firebrick", linetype = "dashed") +
  labs(x = "Pairwise SNP distance", y = "Number of isolate pairs",
       title = "Distances \u2264 100 SNPs (clustering range)") +
  theme_classic(base_size = 11)

if (requireNamespace("gridExtra", quietly = TRUE)) {
  ggsave("snp_distance_histogram.png",
         gridExtra::grid.arrange(p_all, p_zoom, ncol = 2),
         width = 10, height = 4, dpi = 400)
} else {
  ggsave("snp_distance_all.png",  p_all,  width = 5.5, height = 4, dpi = 400)
  ggsave("snp_distance_zoom.png", p_zoom, width = 5.5, height = 4, dpi = 400)
}

sessionInfo()