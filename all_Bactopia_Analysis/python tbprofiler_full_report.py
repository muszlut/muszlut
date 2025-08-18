import json
import sys
import os
import csv

def extract_tbprofiler_full(json_path, threshold=0.9):
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
        dominant = max(lineages, key=lambda x: x.get("fraction", 0.0))
        dom_lineage = dominant.get("lineage", "NA")
        dom_fraction = round(dominant.get("fraction", 0.0), 4)
        dom_family = dominant.get("family", "NA")
        dom_rd = dominant.get("rd", "NA")
        status = "dominant" if dom_fraction >= threshold else "possible_mixed"

        secondary = [
            f"{x.get('lineage','NA')}:{round(x.get('fraction',0.0),4)}"
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

    # -----------------
    # Drug resistance extraction
    # -----------------
    dr_report = data.get("dr_report", {})
    if not dr_report:
        # fallback: check "dr" key if dr_report not present
        dr_report = data.get("dr", {})

    # create dict of drug:status (R=resistant, S=susceptible, NA=unknown)
    resistance = {}
    for drug, info in dr_report.items():
        if isinstance(info, dict):
            resistance[drug] = info.get("prediction", "NA")
        else:
            resistance[drug] = info if info else "NA"

    return sample_id, dom_lineage, dom_fraction, dom_family, dom_rd, status, secondary_lineages, spoligo_octal, spoligo_family, spoligo_sit, resistance

def main(path, output_file=None, mixed_file=None, threshold=0.9):
    all_results = []
    drugs_set = set()

    # Collect JSON files
    json_files = []
    if os.path.isdir(path):
        for root, _, files in os.walk(path):
            for filename in files:
                if filename.endswith(".results.json"):
                    json_files.append(os.path.join(root, filename))
    else:
        json_files.append(path)

    # Extract data
    for f in json_files:
        sample_id, dom_lineage, dom_fraction, dom_family, dom_rd, status, secondary_lineages, spoligo_octal, spoligo_family, spoligo_sit, resistance = extract_tbprofiler_full(f, threshold)
        all_results.append({
            "sample_id": sample_id,
            "dominant_lineage": dom_lineage,
            "fraction": dom_fraction,
            "family": dom_family,
            "RD": dom_rd,
            "status": status,
            "secondary_lineages": secondary_lineages,
            "spoligo_octal": spoligo_octal,
            "spoligo_family": spoligo_family,
            "SIT": spoligo_sit,
            "resistance": resistance
        })
        drugs_set.update(resistance.keys())

    drugs_list = sorted(drugs_set)

    # Prepare header
    header = ["sample_id", "dominant_lineage", "fraction", "family", "RD", "status", "secondary_lineages", "spoligo_octal", "spoligo_family", "SIT"] + drugs_list

    # Write all samples
    if output_file:
        with open(output_file, 'w', newline='') as out:
            writer = csv.writer(out, delimiter='\t')
            writer.writerow(header)
            for r in all_results:
                row = [r["sample_id"], r["dominant_lineage"], r["fraction"], r["family"], r["RD"], r["status"],
                       r["secondary_lineages"], r["spoligo_octal"], r["spoligo_family"], r["SIT"]]
                for drug in drugs_list:
                    row.append(r["resistance"].get(drug, "NA"))
                writer.writerow(row)
    else:
        print("\t".join(header))
        for r in all_results:
            row = [r["sample_id"], r["dominant_lineage"], r["fraction"], r["family"], r["RD"], r["status"],
                   r["secondary_lineages"], r["spoligo_octal"], r["spoligo_family"], r["SIT"]]
            for drug in drugs_list:
                row.append(r["resistance"].get(drug, "NA"))
            print("\t".join(str(c) for c in row))

    # Write mixed samples only
    if mixed_file:
        with open(mixed_file, 'w', newline='') as out:
            writer = csv.writer(out, delimiter='\t')
            writer.writerow(header)
            for r in all_results:
                if r["status"] == "possible_mixed":
                    row = [r["sample_id"], r["dominant_lineage"], r["fraction"], r["family"], r["RD"], r["status"],
                           r["secondary_lineages"], r["spoligo_octal"], r["spoligo_family"], r["SIT"]]
                    for drug in drugs_list:
                        row.append(r["resistance"].get(drug, "NA"))
                    writer.writerow(row)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python tbprofiler_full_report.py <file_or_folder> [all_output.tsv] [mixed_output.tsv] [threshold]")
        sys.exit(1)
    path = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) >= 3 else None
    mixed_file = sys.argv[3] if len(sys.argv) >= 4 else None
    threshold = float(sys.argv[4]) if len(sys.argv) == 5 else 0.9
    main(path, output_file, mixed_file, threshold)
