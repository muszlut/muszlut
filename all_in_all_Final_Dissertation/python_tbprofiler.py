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

# Full standardized drug list you want as columns
drug_list = [
    "rifampicin","isoniazid","ethambutol","pyrazinamide",
    "moxifloxacin","levofloxacin","bedaquiline","delamanid",
    "pretomanid","linezolid","streptomycin","amikacin",
    "kanamycin","capreomycin","clofazimine","ethionamide",
    "para-aminosalicylic_acid","cycloserine"
]

def extract_mutations(variant_list):
    """Return dict: drug -> mutation string"""
    drug_to_mut = {}
    for v in variant_list:
        mut = f"{v['gene']} {v['change']} ({v['freq']})" if "change" in v else v.get("note","")
        for d in v.get("drugs", []):
            drug_to_mut[d["drug"]] = mut
    return drug_to_mut

def main(args):
    # Read sample list or search folder
    if args.samples:
        samples = [x.strip() for x in open(args.samples)]
        files = [os.path.join(args.dir, s, "results", f"{s}.results.json") for s in samples]
    else:
        files = glob.glob(os.path.join(args.dir, "**", f"*{args.suffix}"), recursive=True)
        samples = [os.path.splitext(os.path.basename(f))[0].replace(args.suffix, "") for f in files]

    # Drug sets for resistance-class rules
    FLQ_set = {"moxifloxacin","levofloxacin","ciprofloxacin","ofloxacin"}
    GPA_set = {"bedaquiline","linezolid"}

    # Output CSV
    OUT = open(args.out, "w", newline='')
    fieldnames = [
        "sample","main_lineage","sub_lineage","spoligotype","drtype",
        "target_median_depth","pct_reads_mapped","num_reads_mapped",
        "num_dr_variants","num_other_variants"
    ] + drug_list

    writer = csv.DictWriter(OUT, fieldnames=fieldnames)
    writer.writeheader()

    # Loop samples
    for s, f in tqdm(zip(samples, files), total=len(samples), desc="Processing"):
        data = json.load(open(pp.filecheck(f)))

        # Extract resistance mutations
        dr_mutations = extract_mutations(data.get("dr_variants", []))
        other_mutations = data.get("other_variants", [])

        resistant_drugs = set(dr_mutations.keys())

        # DR-type classification
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

        # Build output row
        row = {
            "sample": s,
            "main_lineage": data.get("main_lineage", "NA"),
            "sub_lineage": data.get("sub_lineage", "NA"),
            "spoligotype": data.get("spoligotype", "NA"),
            "drtype": drtype,
            "target_median_depth": data.get("median_depth", "NA"),
            "pct_reads_mapped": data.get("pct_reads_mapped", "NA"),
            "num_reads_mapped": data.get("num_reads_mapped", "NA"),
            "num_dr_variants": len(data.get("dr_variants", [])),
            "num_other_variants": len(data.get("other_variants", []))
        }

        # Add per-drug mutation (or "-")
        for drug in drug_list:
            row[drug] = dr_mutations.get(drug, "-")

        writer.writerow(row)

    OUT.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extended TBProfiler summary")
    parser.add_argument('--out', required=True, help="Output CSV")
    parser.add_argument('--samples', help="List of sample IDs")
    parser.add_argument('--dir', default="results/", help="Results directory")
    parser.add_argument('--suffix', default=".results.json")
    args = parser.parse_args()
    main(args)
