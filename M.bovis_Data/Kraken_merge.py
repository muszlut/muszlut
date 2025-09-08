import glob
import pandas as pd

# Collect all Kraken2 report files
files = glob.glob("*_kraken2_report.txt")

genus_tables = []
family_tables = []

for f in files:
    # Kraken2 report columns:
    # percent, reads, taxon_reads, rank_code, taxid, name
    df = pd.read_csv(f, sep="\t", header=None, usecols=[0,1,3,4,5],
                     names=["percent", "reads", "rank", "taxid", "name"],
                     engine="python")
    
    sample = f.replace("_kraken2_report.txt", "")
    df["sample"] = sample
    df["name"] = df["name"].str.strip()
    
    # Genus-level
    genus_df = df[df["rank"] == "G"].copy()
    genus_df = genus_df[["name", "reads", "percent", "sample"]]
    genus_tables.append(genus_df)
    
    # Family-level
    family_df = df[df["rank"] == "F"].copy()
    family_df = family_df[["name", "reads", "percent", "sample"]]
    family_tables.append(family_df)

# --- GENUS ---
genus_all = pd.concat(genus_tables)
genus_counts = genus_all.pivot_table(index="name", columns="sample", values="reads", fill_value=0)
genus_perc = genus_all.pivot_table(index="name", columns="sample", values="percent", fill_value=0)
genus_combined = pd.concat({"counts": genus_counts, "percent": genus_perc}, axis=1)
genus_combined.index = ["Genus_" + idx for idx in genus_combined.index]

# --- FAMILY ---
family_all = pd.concat(family_tables)
family_counts = family_all.pivot_table(index="name", columns="sample", values="reads", fill_value=0)
family_perc = family_all.pivot_table(index="name", columns="sample", values="percent", fill_value=0)
family_combined = pd.concat({"counts": family_counts, "percent": family_perc}, axis=1)
family_combined.index = ["Family_" + idx for idx in family_combined.index]

# --- MERGE ---
merged = pd.concat([family_combined, genus_combined])

# --- SORT by mean percent abundance ---
percent_cols = [c for c in merged.columns if isinstance(c, tuple) and c[0] == "percent"]
merged["mean_percent"] = merged[percent_cols].mean(axis=1)
merged = merged.sort_values("mean_percent", ascending=False)

# --- Limit to Top N ---
TOP_N = 30   # adjust as needed
top = merged.head(TOP_N)
others = merged.iloc[TOP_N:]

# Create "Others" row (sum across remaining taxa)
others_counts = others[[c for c in others.columns if isinstance(c, tuple) and c[0] == "counts"]].sum()
others_perc = others[[c for c in others.columns if isinstance(c, tuple) and c[0] == "percent"]].sum()
others_row = pd.concat([others_counts, others_perc])
others_row.name = "Others"

# Final table
merged_top = pd.concat([top.drop(columns="mean_percent"), others_row.to_frame().T])

# Save
merged_top.to_csv("kraken2_top_taxa_family_genus_counts_perc.tsv", sep="\t")

print(f"✅ Done! Created: kraken2_top_taxa_family_genus_counts_perc.tsv (Top {TOP_N} taxa + Others)")
