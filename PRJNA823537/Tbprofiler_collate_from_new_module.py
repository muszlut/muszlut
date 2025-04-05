#! /usr/bin/env python

import sys
import csv
import json
import argparse
import os
from collections import defaultdict
from tqdm import tqdm
import pathogenprofiler as pp
import tbprofiler as tbprofiler

def extract_spoligotype(data):
    """Extract spoligotype from TB-Profiler results JSON."""
    return data.get("spoligotype", "Unknown")

def extract_lineage(data):
    """Extract lineage info from TB-Profiler results JSON."""
    if "lineage" in data and data["lineage"]:
        return data["lineage"].get("lineage", "Unknown")
    return "Unknown"

def extract_amr(data, locus_tag2drugs):
    """Extract AMR mutations and associated genes."""
    amr_list = []
    for var in data.get("dr_variants", []):
        drugs = locus_tag2drugs.get(var["locus_tag"], [])
        amr_list.append((var["gene"], var["change"], ";".join(drugs)))
    return amr_list

def main(args):
    bed_file = f"{sys.base_prefix}/share/tbprofiler/{args.db}.bed"
    locus_tag2drugs = tbprofiler.get_lt2drugs(bed_file)

    if args.samples:
        samples = [x.strip() for x in open(args.samples).readlines()]
    else:
        samples = [x.replace(args.suffix, "") for x in os.listdir(args.dir) if x.endswith(args.suffix)]

    lineage_summary = []
    spoligotype_summary = []
    amr_summary = []

    for s in tqdm(samples):
        filepath = pp.filecheck(f"{args.dir}/{s}{args.suffix}")
        data = json.load(open(filepath))

        lineage = extract_lineage(data)
        spoligotype = extract_spoligotype(data)
        amr_list = extract_amr(data, locus_tag2drugs)

        lineage_summary.append([s, lineage])
        spoligotype_summary.append([s, spoligotype])

        for gene, change, drugs in amr_list:
            amr_summary.append([s, gene, change, drugs])

    # Write lineage summary
    with open(os.path.join(args.dir, "lineage_summary.csv"), 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Sample", "Lineage"])
        writer.writerows(lineage_summary)

    # Write spoligotype summary
    with open(os.path.join(args.dir, "spoligotype_summary.csv"), 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Sample", "Spoligotype"])
        writer.writerows(spoligotype_summary)

    # Write AMR summary
    with open(os.path.join(args.dir, "amr_summary.csv"), 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Sample", "Gene", "Mutation", "Drugs"])
        writer.writerows(amr_summary)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Collate TB-Profiler outputs', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--samples', type=str, help='File with sample names (one per line)')
    parser.add_argument('--dir', default="results/", type=str, help='Directory containing results')
    parser.add_argument('--db', default="tbdb", type=str, help='Database name')
    parser.add_argument('--suffix', default=".results.json", type=str, help='File suffix')
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
