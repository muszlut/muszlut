#!/bin/bash
#SBATCH --job-name=mtb_readlevel
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=1-257%25
#SBATCH --output=/scratch/ma95362/Bio_project_publication/logs/rl.%A_%a.out
#SBATCH --error=/scratch/ma95362/Bio_project_publication/logs/rl.%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu
#
# Read-level presence/absence of the EPTB and T3-ETH candidate genes.
#
# Why: Panaroo calls presence/absence from ASSEMBLIES. group_1905, group_1965,
# group_589, group_1946 and group_2137 are all PPE-family repeats, which
# collapse or fragment during assembly. This maps each isolate's RAW READS to
# the candidate sequences, so presence is measured independently of assembly.
# If a gene is present in the reads of (nearly) all isolates, the Panaroo
# absence calls - and therefore the EPTB association - are assembly artefacts.
#
# Usage:  bash 08a_readlevel_map.sh --check     # build the fastq manifest, map nothing
#         sbatch 08a_readlevel_map.sh           # one array task per isolate
set -euo pipefail

PROJ=/scratch/ma95362/Bio_project_publication
SCRIPTS=${SCRIPTS:-$PROJ/gwas/scripts}
OUT=$PROJ/gwas/readlevel
REF=$OUT/candidates.fa
MAN=$OUT/fastq_manifest.tsv
mkdir -p "$OUT/depth"

module load Mamba/23.11.0-0
source activate mtb_gwas

# ---- 0. candidate sequences (written by step 06) ----
if [ ! -s "$REF" ]; then
    SRC=$PROJ/gwas/followup/followup_sequences.fasta
    [ -s "$SRC" ] || { echo "missing $SRC - run 06_followup.sh first"; exit 1; }
    cp "$SRC" "$REF"
fi

# ---- 1. fastq manifest: sample -> read file(s) ----
if [ ! -s "$MAN" ]; then
    echo "building fastq manifest by searching $PROJ ..."
    find "$PROJ" -name '*.fastq.gz' -o -name '*.fq.gz' 2>/dev/null | grep -v -i 'unpaired\|singleton' > "$OUT/all_fastq.txt" || true
    : > "$MAN"
    while read -r s; do
        files=$(grep -F "$s" "$OUT/all_fastq.txt" | sort | head -2 | tr '\n' ',' | sed 's/,$//')
        [ -n "$files" ] && printf "%s\t%s\n" "$s" "$files" >> "$MAN"
    done < <(tail -n +2 "$SCRIPTS/covariates.tsv" | cut -f1)
    echo "manifest: $(wc -l < "$MAN") of 257 isolates have read files"
fi

if [ "${1:-}" = "--check" ]; then
    echo "--- manifest head"; head -3 "$MAN"
    echo "--- isolates WITHOUT read files:"
    comm -23 <(tail -n +2 "$SCRIPTS/covariates.tsv" | cut -f1 | sort) <(cut -f1 "$MAN" | sort) | head -20
    command -v minimap2 >/dev/null || echo "NOTE: minimap2 missing - run: mamba install -y -n mtb_gwas -c bioconda minimap2"
    echo; echo "if the manifest looks right:  sbatch --array=1-$(wc -l < "$MAN")%25 $SCRIPTS/08a_readlevel_map.sh"
    exit 0
fi

command -v minimap2 >/dev/null || { echo "minimap2 missing: mamba install -y -n mtb_gwas -c bioconda minimap2"; exit 1; }

# ---- 2. one isolate per array task ----
LINE=$(awk -v n="${SLURM_ARRAY_TASK_ID}" 'NR==n' "$MAN")
[ -n "$LINE" ] || { echo "no manifest line ${SLURM_ARRAY_TASK_ID}"; exit 0; }
S=$(echo "$LINE" | cut -f1)
FQ=$(echo "$LINE" | cut -f2 | tr ',' ' ')
echo "$S: $FQ"

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
minimap2 -t 4 -ax sr "$REF" $FQ 2> "$TMP/mm.log" \
  | samtools sort -@ 2 -o "$TMP/$S.bam" -
samtools index "$TMP/$S.bam"

# breadth (fraction of gene covered >= 5x) and mean depth per candidate gene
samtools depth -a -Q 20 "$TMP/$S.bam" \
  | awk -v s="$S" 'BEGIN{OFS="\t"}
      {len[$1]++; sum[$1]+=$3; if($3>=5) cov5[$1]++}
      END{for(g in len) print s, g, len[g], (cov5[g]+0)/len[g], sum[g]/len[g]}' \
  | sort -k2,2 > "$OUT/depth/$S.tsv"

echo "done $S: $(wc -l < "$OUT/depth/$S.tsv") genes"