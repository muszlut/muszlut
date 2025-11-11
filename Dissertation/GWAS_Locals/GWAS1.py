import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# -------------------------------
# Load files
# -------------------------------

# 1. Entropy
entropy = pd.read_csv("alignment_entropy.csv", header=None)
entropy.columns = ["gene", "entropy"]

# 2. Pyseer hits (drug resistance)
pyseer = pd.read_csv("L4.2.2.2_Binary_DR_significant_hits.csv")

# 3. Gene coordinates
coords = pd.read_csv("core-genome.position_cross_reference.txt.gz", sep="\t", header=None)
coords.columns = ["gene", "start", "end"]

# 4. Recombination blocks
recomb = pd.read_csv("core-genome.importation_status.txt", sep="\t")

# -------------------------------
# Merge entropy with coordinates
# -------------------------------
genes = pd.merge(coords, entropy, on="gene", how="left")

# Flag high-entropy genes (top 10%)
entropy_threshold = genes["entropy"].quantile(0.9)
genes["entropy_flag"] = genes["entropy"].apply(lambda x: "High" if x >= entropy_threshold else "Low")

# Merge PySeer hits with coordinates
pyseer_plot = pd.merge(coords, pyseer, left_on="gene", right_on="variant", how="inner")

# -------------------------------
# Plot: Genome overview
# -------------------------------
plt.figure(figsize=(16,4))

# 1. Entropy segments
sns.scatterplot(data=genes, x="start", y="entropy", hue="entropy_flag",
                palette={"High":"red", "Low":"grey"}, s=30)

# 2. Recombination blocks as horizontal blue lines
for idx, row in recomb.iterrows():
    plt.plot([row['Beg'], row['End']], [0.01]*2, color='blue', lw=6, alpha=0.5)

# 3. Pyseer hits as black stars
plt.scatter(pyseer_plot["start"], [0.02]*len(pyseer_plot), marker="*", color="black", s=80, label="GWAS hits")

plt.xlabel("Genome position (bp)")
plt.ylabel("Entropy")
plt.title("High Entropy Genes vs Recombination vs GWAS Hits")
plt.legend(loc="upper right")
plt.tight_layout()
plt.show()

# -------------------------------
# Table of top high-entropy genes
# -------------------------------
top_genes = genes.sort_values(by="entropy", ascending=False).head(30)

# You can manually add gene functions here
top_genes["function"] = ["PDIM biosynthesis" if "pps" in g else
                         "ESX secretion system" if "esp" in g or "mycP" in g else
                         "DNA repair / recombination" if "recC" in g or "ftsK" in g else
                         "Other" for g in top_genes["gene"]]

top_genes.to_csv("Top30_high_entropy_genes.csv", index=False)
print(top_genes[["gene", "entropy", "function"]])
