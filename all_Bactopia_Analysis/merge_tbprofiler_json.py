#!/usr/bin/env python3

import json
import csv
from pathlib import Path

# Set your results directory
RESULTS_DIR = Path("/scratch/ma95362/my_tbprofiler_results")
OUTPUT_FILE = RESULTS_DIR / "merged_tbprofiler_results.csv"

# Prepare CSV header
header = ["ID", "Date", "Strain", "Drug_resistance", "Median_depth", "Lineage", "Spoligotype"]
rows = []

# Loop through sample folders
for sample_dir in RESULTS_DIR.iterdir():
    if not sample_dir.is_dir():
        continue
    
    json_file = sample_dir / "results" / f"{sample_dir.name}.results.json"
    if not json_file.exists():
        continue
    
    with open(json_file) as f:
        data = json.load(f)
    
    # Extract values safely
    ID = data.get("id", "")
    Date = data.get("timestamp", "")
    Strain = data.get("strain", "")
    Drug_resistance = data.get("drug_resistance", "")
    Median_depth = data.get("median_depth", "")

    # Lineage (take the last/highest one in the 'lineage' list if available)
    Lineage = ""
    lineage_list = data.get("lineage", [])
    if lineage_list:
        Lineage = lineage_list[-1].get("name", "")

    # Spoligotype (take first if available)
    Spoligotype = ""
    spoligo_list = data.get("spoligotype", [])
    if spoligo_list:
        Spoligotype = spoligo_list[0].get("pattern", "")

    rows.append([ID, Date, Strain, Drug_resistance, Median_depth, Lineage, Spoligotype])

# Write to CSV
with open(OUTPUT_FILE, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerows(rows)

print(f"Merged results saved to {OUTPUT_FILE}")
# This script merges TBProfiler results from multiple sample directories into a single CSV file.