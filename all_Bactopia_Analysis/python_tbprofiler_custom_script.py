#!/usr/bin/env python

import sys
import csv
import json
import argparse
import os
import glob
from tqdm import tqdm
import pathogenprofiler as pp
import tbprofiler as tbprofiler

def main(args):
    # Path to TBProfiler BED file
    bed_file = os.path.join(sys.base_prefix, "share", "tbprofiler", f"{args.db}.bed")
    locus_tag2drugs = tbprofiler.get_lt2drugs(bed_file)

    # Collect sample files
    if args.samples:
        # If a sample list file is provided
        samples_files = [x.rstrip() for x in open(args.samples).readlines()]
    else:
        # Recursively search for .results.json files in all subdirectories
        samples_files = glob.glob(os.path.join(args.dir, "**", f"*{args.suffix}"), recursive=True)

    if not samples_files:
        print(f"No {args.suffix} files found in {args.dir}")
        return

    # Define drug sets
    FLQ_set = {"moxifloxacin", "levofloxacin", "ciprofloxacin", "ofloxacin"}
    GPA_set = {"bedaquiline", "linezolid"}

    # Open output CSV
    with open(args.out, "w", newline="") as OUT:
        writer = csv.DictWriter(OUT, fieldnames=["sample", "dr-class"])
        writer.writeheader()

        # Process each sample file
        for sample_file in tqdm(samples_files, desc="Processing samples"):
            # Extract sample name from filename
            sample_name = os.path.splitext(os.path.basename(sample_file))[0].replace(args.suffix, "")

            # Load JSON data
            data = json.load(open(pp.filecheck(sample_file)))

            # Collect resistant drugs
            resistant_drugs = {d["drug"] for var in data.get("dr_variants", []) for d in var.get("drugs", [])}

            rif = "rifampicin" in resistant_drugs
            inh = "isoniazid" in resistant_drugs
            flq = bool(FLQ_set.intersection(resistant_drugs))
            gpa = bool(GPA_set.intersection(resistant_drugs))

            # Determine resistance class
            if len(resistant_drugs) == 0:
                drtype = "Sensitive"
            elif (rif and not inh) or (inh and not rif):
                drtype = "Pre-MDR"
            elif (rif and inh) and (not flq and not gpa):
                drtype = "MDR"
            elif (rif and inh) and (flq and not gpa):
                drtype = "Pre-XDR"
            elif (rif and inh) and (flq and gpa):
                drtype = "XDR"
            else:
                drtype = "Other"

            # Write row to CSV
            writer.writerow({"sample": sample_name, "dr-class": drtype})

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TBProfiler merge script", 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--out", type=str, help="Name of CSV output", required=True)
    parser.add_argument("--samples", type=str, help="File with sample paths (optional)")
    parser.add_argument("--dir", default="results/", type=str, help="Directory containing results")
    parser.add_argument("--db", default="tbdb", type=str, help="Database name")
    parser.add_argument("--suffix", default=".results.json", type=str, help="File suffix")
    parser.set_defaults(func=main)

    args = parser.parse_args()
    args.func(args)
