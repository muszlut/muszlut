#!/bin/bash
#SBATCH --job-name=mtb_rl_wg
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=1-257%25
#SBATCH --output=/scratch/ma95362/Bio_project_publication/logs/rlwg.%A_%a.out
#SBATCH --error=/scratch/ma95362/Bio_project_publication/logs/rlwg.%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu
#
# CONTROL for a flaw in steps 08a/08c: those mapped reads against the candidate sequences
# ALONE. With a reference only a few kb long, a read originating from any paralogous PE/PPE
# locus elsewhere in the genome has nowhere else to map, lands on the candidate with high
# MAPQ, and is counted as evidence of presence. The MAPQ filter cannot help, because it only
# resolves ambiguity when the true locus is present in the reference.
#
# Here the reference is the FULL H37Rv genome with the candidate sequences appended as extra
# contigs. Reads from genuine H37Rv paralogues now map to their own locus and compete, so
# depth on a candidate contig reflects reads that match THAT sequence best.
#
# If a candidate's presence call survives this (still ~257/257 for group_2137), the step 08
# conclusion holds. If it collapses, the assembly-artefact claim in the manuscript must be
# rewritten around the whole-genome numbers instead.
#
# Usage:  bash 11_readlevel_wholegenome.sh --prep     # build the reference, map nothing
#         sbatch --array=1-<N>%25 11_readlevel_wholegenome.sh
set -euo pipefail

PROJ=/scratch/ma95362/Bio_project_publication
SCRIPTS=${SCRIPTS:-$PROJ/gwas/scripts}
GBKDIR=/scratch/ma95362/gbk/ncbi_dataset/data/GCF_000195955.2
RL=$PROJ/gwas/readlevel
OUT=$PROJ/gwas/readlevel_wg
REF=$OUT/h37rv_plus_candidates.fa
MAN=$RL/fastq_manifest.tsv
mkdir -p "$OUT/depth"

module load Mamba/23.11.0-0
module load minimap2/2.28-GCCcore-13.2.0
source activate mtb_gwas

# ---- prep: build the combined reference ----
if [ "${1:-}" = "--prep" ] || [ ! -s "$REF" ]; then
    [ -s "$MAN" ] || { echo "missing $MAN - run 08c --check first"; exit 1; }
    CAND=$PROJ/gwas/followup/followup_sequences.fasta
    [ -s "$CAND" ] || { echo "missing $CAND - run 06_followup.sh first"; exit 1; }

    GENOME=$(ls -1 "$GBKDIR"/*.fna 2>/dev/null | head -1 || true)
    if [ -z "$GENOME" ]; then
        echo "no .fna in $GBKDIR - extracting the genome from the GenBank file"
        python3 - "$GBKDIR/genomic.gbk" "$OUT/h37rv.fa" <<'PY'
import re, sys
gbk, out = sys.argv[1], sys.argv[2]
seq, inseq = [], False
with open(gbk) as f:
    for line in f:
        if line.startswith("ORIGIN"): inseq = True; continue
        if inseq:
            if line.startswith("//"): break
            seq.append("".join(re.findall(r"[acgtnACGTN]+", line[10:])))
s = "".join(seq).upper()
if not s: sys.exit("no ORIGIN sequence found")
with open(out, "w") as o:
    o.write(">NC_000962.3 H37Rv\n")
    for i in range(0, len(s), 70): o.write(s[i:i+70] + "\n")
print(f"extracted {len(s)} bp")
PY
        GENOME=$OUT/h37rv.fa
    fi
    echo "genome reference: $GENOME"
    # candidates renamed so they cannot clash with a genome contig name
    awk '/^>/{print ">CAND_" substr($1,2); next} {print}' "$CAND" > "$OUT/candidates_renamed.fa"
    cat "$GENOME" "$OUT/candidates_renamed.fa" > "$REF"
    grep -c '^>' "$REF" | xargs echo "contigs in combined reference:"
    grep '^>' "$OUT/candidates_renamed.fa"
    [ "${1:-}" = "--prep" ] && { echo; echo "now: sbatch --array=1-$(wc -l < "$MAN")%25 $0"; exit 0; }
fi

# ---- one isolate per array task ----
LINE=$(awk -v n="${SLURM_ARRAY_TASK_ID}" 'NR==n' "$MAN")
[ -n "$LINE" ] || { echo "no manifest line ${SLURM_ARRAY_TASK_ID}"; exit 0; }
S=$(echo "$LINE" | cut -f1)
FQ=$(echo "$LINE" | cut -f2 | tr ',' ' ')
echo "$S"

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
minimap2 -t 4 -ax sr "$REF" $FQ 2> "$TMP/mm.log" | samtools sort -@ 2 -o "$TMP/$S.bam" -
samtools index "$TMP/$S.bam"

# depth on the candidate contigs only; MAPQ >= 20 now means "maps here better than anywhere
# else in the genome", which is the claim we actually want to test
samtools depth -a -Q 20 "$TMP/$S.bam" \
  | awk -v s="$S" 'BEGIN{OFS="\t"} $1 ~ /^CAND_/ {
        g=substr($1,6); len[g]++; sum[g]+=$3; if($3>=5) cov5[g]++ }
      END{for(g in len) print s, g, len[g], (cov5[g]+0)/len[g], sum[g]/len[g]}' \
  | sort -k2,2 > "$OUT/depth/$S.tsv"

# genome-wide mean depth, for normalising
gw=$(samtools depth -a -Q 20 -r NC_000962.3 "$TMP/$S.bam" 2>/dev/null \
      | awk '{s+=$3; n++} END{if(n) printf "%.1f", s/n; else print "NA"}')
echo -e "$S\t__genome__\t0\t0\t${gw:-NA}" >> "$OUT/depth/$S.tsv"

echo "done $S: $(wc -l < "$OUT/depth/$S.tsv") rows (genome mean depth ${gw:-NA})"