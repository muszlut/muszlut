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

# Standard columns for drugs
drug_list = [
    "rifampicin","isoniazid","ethambutol","pyrazinamide",
    "moxifloxacin","levofloxacin","bedaquiline","delamanid",
    "pretomanid","linezolid","streptomycin","amikacin",
    "kanamycin","capreomycin","clofazimine","ethionamide",
    "para-aminosalicylic_acid","cycloserine"
]

def extract_mutations(variant_list):
    """Return drug -> mutation string"""
    drug_to_mut = {}
    for v in variant_list:
        mut = f"{v.get('gene','NA')} {v.get('change','NA')} ({v.get('freq','NA')})"
        for d in v.get("drugs", []):
            drug_to_mut[d["drug"]] = mut
    return drug_to_mut

def main(args):

    # ---- Get files ----
    if args.samples:
        samples = [x.strip() for x in open(args.samples)]
        files = [os.path.join(args.dir, s, "results", f"{s}.results.json") for s in samples]
    else:
        files = glob.glob(os.path.join(args.dir, "**", f"*{args.suffix}"), recursive=True)
        samples = [os.path.basename(f).replace(args.suffix, "") for f in files]

    # Drug classification sets
    FLQ_set = {"moxifloxacin","levofloxacin","ciprofloxacin","ofloxacin"}
    GPA_set = {"bedaquiline","linezolid"}

    # ---- Output CSV ----
    OUT = open(args.out, "w", newline='')
    fieldnames = [
        "sample","main_lineage","sub_lineage","spoligotype","drtype",
        "target_median_depth","pct_reads_mapped","num_reads_mapped",
        "num_dr_variants","num_other_variants"
    ] + drug_list
    writer = csv.DictWriter(OUT, fieldnames=fieldnames)
    writer.writeheader()

    # ---- Loop through ALL JSON FILES ----
    for f in tqdm(files, desc="Processing", total=len(files)):

        try:
            data = json.load(open(pp.filecheck(f)))
        except Exception as e:
            print(f"[WARNING] Could not read {f}: {e}")
            continue

        sample = os.path.basename(f).replace(args.suffix, "")

        dr_mutations = extract_mutations(data.get("dr_variants", []))
        resistant_drugs = set(dr_mutations.keys())

        # ---- DR classification ----
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

        # ---- Prepare row ----
        row = {
            "sample": sample,
            "main_lineage": data.get("main_lineage", "NA"),
            "sub_lineage": data.get("sub_lineage", "NA"),
            "spoligotype": data.get("spoligotype", "NA"),
            "drtype": drtype,
            "target_median_depth": data.get("median_depth", "NA"),
            "pct_reads_mapped": data.get("pct_reads_mapped", "NA"),
            "num_reads_mapped": data.get("num_reads_mapped", "NA"),
            "num_dr_variants": len(data.get("dr_variants", [])),
            "num_other_variants": len(data.get("other_variants", [])),
        }

        for drug in drug_list:
            row[drug] = dr_mutations.get(drug, "-")

        writer.writerow(row)

    OUT.close()
    print(f"\n[DONE] Wrote output file: {args.out}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extended TBProfiler summary")
    parser.add_argument('--out', required=True, help="Output CSV")
    parser.add_argument('--samples', help="List of sample IDs")
    parser.add_argument('--dir', default="results/", help="Results directory")
    parser.add_argument('--suffix', default=".results.json")
    args = parser.parse_args()
    main(args)
