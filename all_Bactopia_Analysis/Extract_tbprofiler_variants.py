#!/usr/bin/env python3
import os
import sys
import json
import glob
import pandas as pd

def parse_tbprofiler_results(json_file):
    with open(json_file) as f:
        data = json.load(f)

    sample_id = os.path.basename(json_file).replace(".results.json", "")
    records = []

    # Drug resistance mutations
    if "dr_variants" in data:
        for var in data["dr_variants"]:
            record = {
                "Sample": sample_id,
                "Gene": var.get("gene", ""),
                "Change": var.get("change", ""),
                "Type": var.get("type", ""),
                "Nucleotide Change": var.get("nucleotide_change", ""),
                "Amino Acid Change": var.get("protein_change", ""),
                "Drugs": ",".join(var.get("drug", [])) if isinstance(var.get("drug", []), list) else var.get("drug", ""),
                "Confidence": var.get("confidence", ""),
                "Freq": var.get("freq", ""),
                "Depth": var.get("depth", ""),
            }
            records.append(record)

    # Other variants (non-drug resistance)
    if "other_variants" in data:
        for var in data["other_variants"]:
            record = {
                "Sample": sample_id,
                "Gene": var.get("gene", ""),
                "Change": var.get("change", ""),
                "Type": var.get("type", ""),
                "Nucleotide Change": var.get("nucleotide_change", ""),
                "Amino Acid Change": var.get("protein_change", ""),
                "Drugs": "None",
                "Confidence": "NA",
                "Freq": var.get("freq", ""),
                "Depth": var.get("depth", ""),
            }
            records.append(record)

    return records


def main(results_dir, output_tsv):
    all_records = []

    # Search recursively for *.results.json
    for json_file in glob.glob(os.path.join(results_dir, "**", "*.results.json"), recursive=True):
        all_records.extend(parse_tbprofiler_results(json_file))

    if not all_records:
        print("No variants found in results.")
        return

    df = pd.DataFrame(all_records)
    df.to_csv(output_tsv, sep="\t", index=False)
    print(f"[INFO] Extracted {len(all_records)} variant records into {output_tsv}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python Extract_variants.py <results_dir> <output.tsv>")
        sys.exit(1)

    results_dir = sys.argv[1]
    output_tsv = sys.argv[2]
    main(results_dir, output_tsv)
