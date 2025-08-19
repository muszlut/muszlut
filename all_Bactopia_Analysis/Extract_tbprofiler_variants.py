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
    for var in data.get("dr_variants", []):
        gene_name = var.get("gene") or var.get("gene_name") or var.get("locus_tag") or "NA"
        aa_change = var.get("protein_change") or var.get("change") or "NA"
        nuc_change = var.get("nucleotide_change") or "NA"

        drugs = var.get("drug") or var.get("drugs") or []
        # Handle list of dicts, strings, or None values
        if isinstance(drugs, list):
            drugs_str = ",".join([d.get("name") if isinstance(d, dict) and d else str(d) for d in drugs if d])
        else:
            drugs_str = str(drugs)

        record = {
            "Sample": sample_id,
            "Gene": gene_name,
            "Change": aa_change,
            "Type": var.get("type", ""),
            "Nucleotide Change": nuc_change,
            "Amino Acid Change": aa_change,
            "Drugs": drugs_str,
            "Confidence": var.get("confidence") or var.get("conf") or "",
            "Freq": var.get("freq") or var.get("estimated_fraction") or "",
            "Depth": var.get("depth", "")
        }
        records.append(record)

    # Other variants (non-drug resistance)
    for var in data.get("other_variants", []):
        gene_name = var.get("gene") or var.get("gene_name") or var.get("locus_tag") or "NA"
        aa_change = var.get("protein_change") or var.get("change") or "NA"
        nuc_change = var.get("nucleotide_change") or "NA"

        record = {
            "Sample": sample_id,
            "Gene": gene_name,
            "Change": aa_change,
            "Type": var.get("type", ""),
            "Nucleotide Change": nuc_change,
            "Amino Acid Change": aa_change,
            "Drugs": "None",
            "Confidence": "NA",
            "Freq": var.get("freq") or var.get("estimated_fraction") or "",
            "Depth": var.get("depth", "")
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
        print("Usage: python Extract_variants_fixed3.py <results_dir> <output.tsv>")
        sys.exit(1)

    results_dir = sys.argv[1]
    output_tsv = sys.argv[2]
    main(results_dir, output_tsv)
