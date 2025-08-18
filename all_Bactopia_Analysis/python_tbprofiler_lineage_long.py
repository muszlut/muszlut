import json
import sys
import os
import csv

def extract_lineage_fractions(json_path):
with open(json_path, 'r') as f:
data = json.load(f)
sample_id = data.get("id", os.path.basename(json_path).replace(".results.json", ""))

output = []
for lineage_data in data.get("lineage", []):
lineage = lineage_data.get("lineage", "NA")
fraction = round(lineage_data.get("fraction", 0.0), 4)
output.append((sample_id, lineage, fraction))
return output

def main(path, output_file=None):
results = []

if os.path.isdir(path):
for root, _, files in os.walk(path):
for filename in files:
if filename.endswith(".results.json"):
full_path = os.path.join(root, filename)
results.extend(extract_lineage_fractions(full_path))
else:
results.extend(extract_lineage_fractions(path))

# Write output
if output_file:
with open(output_file, 'w', newline='') as out:
writer = csv.writer(out, delimiter='\t')
writer.writerow(["sample_id", "lineage", "fraction"])
writer.writerows(results)
else:
print("sample_id\tlineage\tfraction")
for row in results:
print("\t".join(str(col) for col in row))

if __name__ == "__main__":
if len(sys.argv) < 2:
print("Usage: python tbprofiler_lineage_long.py <file_or_folder> [output.tsv]")
sys.exit(1)
path = sys.argv[1]
output_file = sys.argv[2] if len(sys.argv) == 3 else None
main(path, output_file)