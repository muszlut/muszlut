#!/bin/bash
#SBATCH --job-name=mtb_rl_dualmapq
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=1-257%25
#SBATCH --output=/scratch/ma95362/Bio_project_publication/logs/rlq.%A_%a.out
#SBATCH --error=/scratch/ma95362/Bio_project_publication/logs/rlq.%A_%a.err
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
# Usage:  bash 08c_readlevel_map_dualmapq.sh --check     # build the fastq manifest, map nothing
#         sbatch 08c_readlevel_map_dualmapq.sh           # one array task per isolate
set -euo pipefail

PROJ=/scratch/ma95362/Bio_project_publication
SCRIPTS=${SCRIPTS:-$PROJ/gwas/scripts}
OUT=$PROJ/gwas/readlevel
REF=$OUT/candidates.fa
MAN=$OUT/fastq_manifest.tsv
mkdir -p "$OUT/depth"

module load Mamba/23.11.0-0
# minimap2 comes from the cluster module, not the conda env (loaded BEFORE activating so
# that conda's samtools still takes precedence on PATH).
module load minimap2/2.28-GCCcore-13.2.0
source activate mtb_gwas

# ---- 0. candidate sequences (written by step 06) ----
if [ ! -s "$REF" ]; then
    SRC=$PROJ/gwas/followup/followup_sequences.fasta
    [ -s "$SRC" ] || { echo "missing $SRC - run 06_followup.sh first"; exit 1; }
    cp "$SRC" "$REF"
fi

# ---- 1. fastq manifest: sample -> read file(s) ----
# Built ONLY by --check: array tasks must never build it concurrently (they would
# all append to the same file). Array tasks read it and fail if it is missing.
if [ "${1:-}" != "--check" ] && [ ! -s "$MAN" ]; then
    echo "missing $MAN - run 'bash $0 --check' first"; exit 1
fi
if [ "${1:-}" = "--check" ]; then
    echo "building fastq manifest by searching $PROJ ..."
    rm -f "$MAN"
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
    echo; echo "if the manifest looks right:  sbatch --array=1-$(wc -l < "$MAN")%25 $0"
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

# breadth (fraction of gene covered >= 5x) and mean depth per candidate gene.
#
# Two versions, because MAPQ matters for repeats: these are PPE genes, so reads map
# ambiguously and get MAPQ 0. Filtering at -Q 20 therefore DISCARDS true reads and can
# make a present gene look absent. -Q 0 keeps them but cannot tell a true locus from a
# paralogue. Neither is right on its own; a gene present in both is solidly present, and
# a gene whose call flips between them is unreliable - which is itself the finding.
mkdir -p "$OUT/depth" "$OUT/depth_q0"
for q in 20 0; do
    d=$([ "$q" = 20 ] && echo "$OUT/depth" || echo "$OUT/depth_q0")
    samtools depth -a -Q "$q" "$TMP/$S.bam" \
      | awk -v s="$S" 'BEGIN{OFS="\t"}
          {len[$1]++; sum[$1]+=$3; if($3>=5) cov5[$1]++}
          END{for(g in len) print s, g, len[g], (cov5[g]+0)/len[g], sum[g]/len[g]}' \
      | sort -k2,2 > "$d/$S.tsv"
done

echo "done $S: $(wc -l < "$OUT/depth/$S.tsv") genes (Q20), $(wc -l < "$OUT/depth_q0/$S.tsv") (Q0)"