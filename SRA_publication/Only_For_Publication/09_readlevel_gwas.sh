#!/bin/bash
#SBATCH --job-name=mtb_rl_gwas
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/ma95362/Bio_project_publication/logs/log.%j.out
#SBATCH --error=/scratch/ma95362/Bio_project_publication/logs/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu
#
# The like-for-like test: run the SAME pyseer LMM used in step 03, but with the candidate
# genes genotyped from READS instead of from assemblies.
#
# Step 08's Fisher tests are not comparable to the pyseer p-values because Fisher ignores
# population structure. This closes that gap: identical kinship matrix, identical
# phenotype, identical model - only the genotype source differs. If the assembly-based
# p-value collapses when reads are used, the association was an assembly artefact.
#
# Also tests the T3-ETH phenotypes, to confirm group_2137 is not a clade marker.
set -euo pipefail

module load Mamba/23.11.0-0
source activate mtb_gwas

PROJ=/scratch/ma95362/Bio_project_publication
SCRIPTS=${SCRIPTS:-$PROJ/gwas/scripts}
PY=$PROJ/gwas/pyseer
RL=$PROJ/gwas/readlevel
OUT=$RL/gwas
mkdir -p "$OUT"; cd "$OUT"

for f in "$SCRIPTS/make_readlevel_rtab.py" "$SCRIPTS/phenotypes.tsv" "$SCRIPTS/covariates.tsv" \
         "$PY/phylogeny_K.tsv"; do
    [ -s "$f" ] || { echo "missing input: $f"; exit 1; }
done
[ -d "$RL/depth" ] || { echo "missing $RL/depth - run 08a first"; exit 1; }

echo "===== 1. read-level genotype matrix"
python3 "$SCRIPTS/make_readlevel_rtab.py" "$RL/depth" "$SCRIPTS/covariates.tsv" "$OUT/readlevel"

npat () { sort -u "$1" | wc -l; }

echo; echo "===== 2. pyseer LMM on read-level genotypes"
# min-af 0 / max-af 1: only 5 genes here, and genes present in all isolates are
# informative (they are the proof of artefact), so nothing is filtered out.
run () {   # $1 phenotype column  $2 phenotype file
    local tag="$1_readlevel"
    pyseer --phenotypes "$2" --phenotype-column "$1" --pres "$OUT/readlevel.Rtab" \
           --similarity "$PY/phylogeny_K.tsv" --lmm \
           --min-af 0 --max-af 1 --cpu 4 \
           --output-patterns "patterns_${tag}.txt" > "pyseer_${tag}.tsv" 2> "pyseer_${tag}.log"
    local n thr
    n=$(npat "patterns_${tag}.txt"); thr=$(awk -v n="$n" 'BEGIN{printf "%.3g", 0.05/n}')
    echo "--- $1 (patterns=$n, threshold=$thr)"
    column -t -s$'\t' "pyseer_${tag}.tsv" | cut -c1-150
    echo
}
run EPTB "$SCRIPTS/phenotypes.tsv"
for c in T3ETH_vs_otherL4222 T3ETH_vs_otherL4 T3ETH_vs_all; do
    [ -s "$SCRIPTS/pheno_$c.tsv" ] && run "$c" "$SCRIPTS/pheno_$c.tsv"
done

echo "===== 3. side by side with the assembly-based results"
printf "%-34s %-14s %-14s\n" variant p_assembly p_readlevel
for v in group_1905 group_1965 group_589 group_1946 group_2137; do
    pa=$(awk -F'\t' -v v="$v" '$1==v{print $4}' "$PY/pyseer_EPTB_genes.tsv" 2>/dev/null | head -1)
    pr=$(awk -F'\t' -v v="$v" '$1==v{print $4}' "pyseer_EPTB_readlevel.tsv" 2>/dev/null | head -1)
    printf "%-34s %-14s %-14s\n" "$v" "${pa:-not tested}" "${pr:-NA}"
done

echo; echo "===== done ($(date)); results: $OUT"