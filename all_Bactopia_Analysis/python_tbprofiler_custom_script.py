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
    # Prepare locus tag -> drugs mapping
    bed_file = os.path.join(sys.base_prefix, "share", "tbprofiler", f"{args.db}.bed")
    locus_tag2drugs = tbprofiler.get_lt2drugs(bed_file)

    # Get samples: either from file or by scanning the directory recursively
    if args.samples:
        samples = [x.strip() for x in open(args.samples)]
        files = [os.path.join(args.dir, s, "results", f"{s}.results.json") for s in samples]
    else:
        files = glob.glob(os.path.join(args.dir, "**", f"*{args.suffix}"), recursive=True)
        samples = [os.path.splitext(os.path.basename(f))[0].replace(args.suffix, "") for f in files]

    FLQ_set = {"moxifloxacin","levofloxacin","ciprofloxacin","ofloxacin"}
    GPA_set = {"bedaquiline","linezolid"}

    # Open output CSV
    OUT = open(args.out, "w", newline='')
    writer = csv.DictWriter(
        OUT,
        fieldnames=[
            "ID","Date","Strain","Drug_resistance","dr-class",
            "Median_depth","Lineage","Sublineage","Spoligotype"
        ]
    )
    writer.writeheader()

    # Process each sample
    for s, f in tqdm(zip(samples, files), total=len(samples), desc="Processing samples"):
        data = json.load(open(pp.filecheck(f)))
        resistant_drugs = set()
        for var in data.get("dr_variants", []):
            for d in var.get("drugs", []):
                resistant_drugs.add(d["drug"])

        rif = "rifampicin" in resistant_drugs
        inh = "isoniazid" in resistant_drugs
        flq = len(FLQ_set.intersection(resistant_drugs)) > 0
        gpa = len(GPA_set.intersection(resistant_drugs)) > 0

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

        writer.writerow({
            "ID": data.get("id", s),
            "Date": data.get("date", "NA"),
            "Strain": data.get("strain", "NA"),
            "Drug_resistance": ";".join(sorted(resistant_drugs)) if resistant_drugs else "None",
            "dr-class": drtype,
            "Median_depth": data.get("median_depth", "NA"),
            "Lineage": data.get("main_lineage", "NA"),
            "Sublineage": data.get("sub_lineage", "NA"),
            "Spoligotype": data.get("spoligotype", "NA")
        })

    OUT.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='TBProfiler extended summary script', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--out', type=str, required=True, help='Name of CSV output')
    parser.add_argument('--samples', type=str, help='File listing samples (optional)')
    parser.add_argument('--dir', default="results/", type=str, help='Directory containing results')
    parser.add_argument('--db', default="tbdb", type=str, help='Database name')
    parser.add_argument('--suffix', default=".results.json", type=str, help='File suffix')
    args = parser.parse_args()
    main(args)
