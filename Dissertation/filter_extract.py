import pandas as pd
 
# Load your GWAS results
df = pd.read_csv("T3ETH_gwas.txt", sep="\t")

# Remove "bad-chisq" entries if they exist
df = df[~df["notes"].astype(str).str.contains("bad-chisq", na=False)]

# Define Bonferroni-corrected threshold
sig_threshold = 1.21e-4

# Filter for significant hits
sig_hits = df[df["lrt-pvalue"] < sig_threshold]

# Save significant hits
sig_hits.to_csv("T3-ETH_significant_hits.csv", index=False)

print(f"Total significant hits: {len(sig_hits)}")
print(f"Bonferroni threshold: {sig_threshold}")
print("Significant hits saved to T3-ETH_significant_hits.csv")
print(sig_hits)
import os