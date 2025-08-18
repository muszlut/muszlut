import json
import sys
import os
import csv
import statistics

def extract_dr_report(json_path):
    with open(json_path, 'r') as f:
        data = json.load(f)
    sample_id = data.get("id", os.path.basename(json_path).replace(".results.json", ""))

    # Resistance summary
    dr_report = data.get("dr_report", {})
    if not dr_report:
        dr_report = data.get("dr", {})

    # Variant-level info
    variants = data.get("variants", [])

    output_rows = []

    for drug, info in dr_report.items():
        if isinstance(info, dict):
            status = info.get("prediction", "NA")
        else:
            status = info if info else "NA"

        # Filter variants related to this drug
        drug_variants = [v for v in variants if drug in v.get("gene_associated_drug", [])]

        # Calculate median depth of resistant variants (if available)
        depths = [v.get("depth", 0) for v in drug_variants if v.get("depth") is not None]
        median_depth = round(statistics.median(depths), 2) if depths else "NA"

        if drug_variants:
            for v in drug_variants:
                row = [
                    sample_id,
                    drug,
                    status,
                    median_depth,
                    v.get("genome_position", "NA"),
                    v.get("locus_tag", "NA"),
                    v.get("gene_name", "NA"),
                    v.get("variant_type", "NA"),
                    v.get("change", "NA"),
                    v.get("depth", "NA"),
                    v.get("estimated_fraction", "NA"),
                    ",".join(v.get("gene_associated_drug", [])),
                    v.get("confidence", "NA"),
                    v.get("comment", "NA")
                ]
                output_rows.append(row)
        else:
            # No variants, still include summary row
            row = [
                sample_id,
                drug,
                status,
                median_depth,
                "NA","NA","NA","NA","NA","NA","NA",drug,"NA","NA"
            ]
            output_rows.append(row)

    return output_rows

def main(path, output_file="tbprofiler_dr_report.tsv"):
    all_rows = []

    # Collect JSON files
    json_files = []
    if os.path.isdir(path):
        for root, _, files in os.walk(path):
            for filename in files:
                if filename.endswith(".results.json"):
                    json_files.append(os.path.join(root, filename))
    else:
        json_files.append(path)

    for f in json_files:
        all_rows.extend(extract_dr_report(f))

    # Header
    header = ["Sample_ID", "Drug", "Resistance", "Median_Depth_Resistant_Variants",
              "Genome_Position", "Locus_Tag", "Gene_Name", "Variant_Type", "Change",
              "Depth", "Estimated_Fraction", "Gene_Associated_Drug", "Confidence", "Comment"]

    # Write TSV
    with open(output_file, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(header)
        writer.writerows(all_rows)

    print(f"Drug resistance report saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python tbprofiler_drug_report.py <file_or_folder> [output_file.tsv]")
        sys.exit(1)
    path = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) >=3 else "tbprofiler_dr_report.tsv"
    main(path, output_file)
