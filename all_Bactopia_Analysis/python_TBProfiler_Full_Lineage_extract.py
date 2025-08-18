import json
import sys
import os
import csv

def extract_tbprofiler_info(json_path, threshold=0.9):
    with open(json_path, 'r') as f:
        data = json.load(f)
    sample_id = data.get("id", os.path.basename(json_path).replace(".results.json", ""))

    # -----------------
    # Lineage extraction
    # -----------------
    lineages = data.get("lineage", [])
    if not lineages:
        dom_lineage = "NA"
        dom_fraction = 0.0
        dom_family = "NA"
        dom_rd = "NA"
        status = "no_lineage"
        secondary_lineages = "NA"
    else:
        # Dominant lineage
        dominant = max(lineages, key=lambda x: x.get("fraction", 0.0))
        dom_lineage = dominant.get("lineage", "NA")
        dom_fraction = round(dominant.get("fraction", 0.0), 4)
        dom_family = dominant.get("family", "NA")
        dom_rd = dominant.get("rd", "NA")

        # Status
        status = "dominant" if dom_fraction >= threshold else "possible_mixed"

        # Secondary lineages (lineage: fraction)
        secondary = [
            f"{x.get('lineage', 'NA')}:{round(x.get('fraction',0.0),4)}"
            for x in lineages if x.get("lineage") != dom_lineage and x.get("fraction",0.0) > 0
        ]
        secondary_lineages = ";".join(secondary) if secondary else "NA"

    # -----------------
    # Spoligotype extraction
    # -----------------
    spoligo_data = data.get("spoligotype", {})
    spoligo_octal = spoligo_data.get("octal", "NA")
    spoligo_family = spoligo_data.get("family", "NA")
    spoligo_sit = spoligo_data.get("SIT", "NA")

    return [(sample_id, dom_lineage, dom_fraction, dom_family, dom_rd, status, 
             secondary_lineages, spoligo_octal, spoligo_family, spoligo_sit)]

def main(path, output_file=None, mixed_file=None, threshold=0.9):
    results = []

    # Collect JSON files
    if os.path.isdir(path):
        for root, _, files in os.walk(path):
            for filename in files:
                if filename.endswith(".results.json"):
                    full_path = os.path.join(root, filename)
                    results.extend(extract_tbprofiler_info(full_path, threshold))
    else:
        results.extend(extract_tbprofiler_info(path, threshold))

    # -----------------
    # Output headers
    # -----------------
    header = ["sample_id", "dominant_lineage", "fraction", "family", "RD", 
              "status", "secondary_lineages", "spoligo_octal", "spoligo_family", "SIT"]

    # Write all samples
    if output_file:
        with open(output_file, 'w', newline='') as out:
            writer = csv.writer(out, delimiter='\t')
            writer.writerow(header)
            writer.writerows(results)
    else:
        print("\t".join(header))
        for row in results:
            print("\t".join(str(col) for col in row))

    # Write only mixed samples
    if mixed_file:
        mixed_results = [r for r in results if r[5] == "possible_mixed"]
        with open(mixed_file, 'w', newline='') as out:
            writer = csv.writer(out, delimiter='\t')
            writer.writerow(header)
            writer.writerows(mixed_results)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python tbprofiler_full_extract.py <file_or_folder> [all_output.tsv] [mixed_output.tsv] [threshold]")
        sys.exit(1)
    path = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) >= 3 else None
    mixed_file = sys.argv[3] if len(sys.argv) >= 4 else None
    threshold = float(sys.argv[4]) if len(sys.argv) == 5 else 0.9
    main(path, output_file, mixed_file, threshold)
