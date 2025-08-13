#!/usr/bin/env python3

import pandas as pd
import json
import glob
import os

# Set your TBProfiler results directory
OUTDIR = "/scratch/ma95362/my_tbprofiler_results"

# Find all JSON result files in all sample subdirectories
json_files = glob.glob(os.path.join(OUTDIR, "*/results/*.results.json"))

# Prepare a list to store each sample's info
data = []

for f in json_files:
    with open(f) as jfile:
        report = json.load(jfile)
        
        # Extract fields you want (customize as needed)
        data.append({
            "ID": report.get("ID"),
            "Date": report.get("Date"),
            "Strain": report.get("Strain"),
            "Drug_resistance": report.get("Drug-resistance"),
            "Median_depth": report.get("Median Depth"),
            # Example: extract top-level lineage summary
            "Lineage": report.get("Lineage", {}).get("lineage", ""),
            # Optionally, include spoligotype if available
            "Spoligotype": report.get("Spoligotype", {}).get("octal", "")
        })

# Convert list of dicts to a pandas DataFrame
df = pd.DataFrame(data)

# Save merged results to CSV
merged_csv = os.path.join(OUTDIR, "merged_tbprofiler_results.csv")
df.to_csv(merged_csv, index=False)

print(f"Merged results saved to {merged_csv}")
print(df.head())
