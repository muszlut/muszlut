#!/bin/bash
#SBATCH --job-name=pyseer_gwas
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Load Python if not automatically available
module load python
PANAROO_DIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_of_1368/panaroo"
cd "${PANAROO_DIR}" || { echo "❌ ERROR: Panaroo directory not found."; exit 1; }
# === Define file paths ===
PPA="/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_of_1368/panaroo/gene_presence_absence.csv"
META="/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_of_1368/Full_metadata.tab"
OUT1="cysE_presence_absence.tsv"
OUT2="cysE_presence_binary.tsv"

echo "Extracting cysE-group_2002-group_2003 presence/absence from: $PPA"

# === Step 1: Extract header and target row ===
head -n1 "$PPA" > header.tmp
grep -P "^cysE-group_2002-group_2003\b" "$PPA" > "$OUT1"

# Fallback: use a loose search if exact name not found
if [ ! -s "$OUT1" ]; then
  echo "Exact name not found, using loose search for cysE..."
  grep -i "cysE" "$PPA" > "$OUT1"
fi

echo "Partial results preview:"
cut -f1-10 "$OUT1" | column -t -s$'\t'

# === Step 2: Generate binary (1/0) presence/absence table ===
python3 - <<'PY'
import csv, sys
pfile = "panaroo_output/gene_presence_absence.csv"
meta = "metadata.tsv"  # isolate\tlineage
out = "cysE_presence_binary.tsv"
target = "cysE-group_2002-group_2003"

# Load header to get isolate names
with open(pfile) as f:
    cols = next(f).strip().split(',')
    isolate_names = cols[14:] if len(cols) > 14 else cols[6:]  # adjust if needed

# Find the target row
with open(pfile) as f:
    for line in f:
        if target in line.split(',')[0]:
            row = line.strip().split(',')
            break
    else:
        print("Target not found; showing available lines containing cysE:")
        with open(pfile) as f2:
            for l in f2:
                if "cysE" in l:
                    print(l)
        sys.exit(1)

# Mark presence/absence (1 = present, 0 = absent)
presence = []
for cell in row[14:]:
    presence.append("1" if cell.strip() else "0")

# Write results
with open(out, "w") as fo:
    fo.write("isolate\tpresence\n")
    for iso, p in zip(isolate_names, presence):
        fo.write(f"{iso}\t{p}\n")
print("✅ Wrote", out)
PY

# === Step 3: Join with lineage metadata and summarize counts ===
echo "Joining metadata with presence table..."
join -t $'\t' -1 1 -2 1 <(sort -k1,1 "$META") <(sort -k1,1 "$OUT2") | \
  awk -F '\t' '{counts[$2","$3]++} END{for(k in counts) print k,counts[k]}' | sort > lineage_cysE_counts.tsv

echo "✅ Summary written to lineage_cysE_counts.tsv"
echo "Done."
