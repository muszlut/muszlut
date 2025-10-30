import pandas as pd
import os
 
# Load your GWAS results
df = pd.read_csv("EPTB_vs_PTb_gwas.txt", sep="\t")

# Remove "bad-chisq" entries if they exist
df = df[~df["notes"].astype(str).str.contains("bad-chisq", na=False)]

# Define Bonferroni-corrected threshold
sig_threshold = 1.06E-04

# Filter for significant hits
sig_hits = df[df["lrt-pvalue"] < sig_threshold]

# Save significant hits
sig_hits.to_csv("EPTB_vs_PTb_significant_hits.csv", index=False)

print(f"Total significant hits: {len(sig_hits)}")
print(f"Bonferroni threshold: {sig_threshold}")
print("Significant hits saved to EPTB_vs_PTb_significant_hits.csv")
print(sig_hits)
#End