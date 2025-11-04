#!/bin/bash
#SBATCH --job-name=cysE_presence_analysis
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/cysE_presence_%j.out
#SBATCH --error=/scratch/ma95362/scratch/cysE_presence_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# === Load Python ===
module load Python/3.13.1-GCCcore-14.2.0

# === Define paths ===
PANAROO_DIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_of_1368/panaroo"
PPA="$PANAROO_DIR/gene_presence_absence.csv"
META="/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_of_1368/Full_metadata.tab"
OUT1="$PANAROO_DIR/cysE_presence_absence.tsv"
OUT2="$PANAROO_DIR/cysE_presence_binary.tsv"

cd "$PANAROO_DIR" || { echo "❌ ERROR: Panaroo directory not found."; exit 1; }

echo "🔍 Extracting cysE-group_2002-group_2003 from: $PPA"

# Step 1: Extract the header and matching row
head -n1 "$PPA" > header.tmp
grep -P "^cysE-group_2002-group_2003\b" "$PPA" > "$OUT1"

# Fallback: loose search
if [ ! -s "$OUT1" ]; then
  echo "⚠️ Exact match not found; searching loosely for 'cysE'..."
  grep -i "cysE" "$PPA" > "$OUT1"
fi

# Preview
echo "🧾 Partial output preview:"
cut -f1-10 "$OUT1" | column -t -s$'\t'

# Step 2: Convert to binary presence/absence table
python3 - <<PY
import sys, csv

pfile = "$PPA"
out = "$OUT2"
target = "cysE-group_2002-group_2003"

# Load header
with open(pfile) as f:
    first_line = f.readline()
    delimiter = ',' if ',' in first_line else '\t'
    cols = first_line.strip().split(delimiter)
    isolate_names = cols[14:] if len(cols) > 14 else cols[6:]

# Find target line
row = None
with open(pfile) as f:
    for line in f:
        if target in line:
            row = line.strip().split(delimiter)
            break

if row is None:
    print("❌ Target not found; available lines with 'cysE':")
    with open(pfile) as f2:
        for l in f2:
            if "cysE" in l:
                print(l)
    sys.exit(1)

# Create binary list
presence = ["1" if cell.strip() else "0" for cell in row[14:]]

# Write output
with open(out, "w") as fo:
    fo.write("isolate\tpresence\n")
    for iso, p in zip(isolate_names, presence):
        fo.write(f"{iso}\t{p}\n")

print(f"✅ Wrote binary presence table: {out}")
PY

# Step 3: Merge with metadata
echo "🔗 Joining metadata and summarizing..."
join -t $'\t' -1 1 -2 1 <(sort -k1,1 "$META") <(sort -k1,1 "$OUT2") | \
  awk -F '\t' '{counts[$2","$3]++} END{for(k in counts) print k,counts[k]}' | sort > lineage_cysE_counts.tsv

echo "✅ Done! Summary saved to lineage_cysE_counts.tsv"
