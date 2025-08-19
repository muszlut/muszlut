#!/usr/bin/env python3
import os
import sys
import json
import glob

def parse_tbprofiler_json(json_path):
    with open(json_path, "r") as f:
        data = json.load(f)
    
    sample = data.get("sample", os.path.basename(json_path).replace(".results.json", ""))
    median_depth = data.get("median_depth", "NA")

    rows = []

    # Drug resistance variants
    for var in data.get("dr_variants", []):
        rows.append([
            sample,
            var.get("genome_pos", "NA"),
            var.get("locus_tag", "NA"),
            var.get("gene", "NA"),
            var.get("var_type", "NA"),
            var.get("change", "NA"),
            var.get("depth", "NA"),
            var.get("freq", "NA"),
            ",".join(var.get("drugs", [])),
            var.get("conf", "NA"),
            var.get("comment", ""),
            median_depth
        ])

    # Other variants
    for var in data.get("other_variants", []):
        rows.append([
            sample,
            var.get("genome_pos", "NA"),
            var.get("locus_tag", "NA"),
            var.get("gene", "NA"),
            var.get("var_type", "NA"),
            var.get("change", "NA"),
            var.get("depth", "NA"),
            var.get("freq", "NA"),
            "",  # No drug association
            "",  # No confidence
            "",  # No comment
            median_depth
        ])

    return rows


def main(results_dir, output_file):
    all_rows = []
    header = [
        "Sample", "Genome Position", "Locus Tag", "Gene Name", "Variant Type", 
        "Change", "Depth", "Estimated Fraction", "Gene Associated Drug", 
        "Confidence", "Comment", "Median Depth"
    ]

    json_files = glob.glob(os.path.join(results_dir, "*", "*.results.json"))
    for jf in json_files:
        all_rows.extend(parse_tbprofiler_json(jf))

    with open(output_file, "w") as out:
        out.write("\t".join(header) + "\n")
        for row in all_rows:
            out.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python extract_tbprofiler_variants.py <results_dir> <output.tsv>")
        sys.exit(1)

    results_dir = sys.argv[1]
    output_file = sys.argv[2]
    main(results_dir, output_file)
