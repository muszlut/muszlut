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

# Which depth set to genotype from, by argument or environment variable:
#   depth     MAPQ >= 20 (default)      -> results in gwas/readlevel/gwas
#   depth_q0  unfiltered, from 08c      -> results in gwas/readlevel/gwas_q0
DSET=${1:-${DSET:-depth}}
SUF=$([ "$DSET" = "depth" ] && echo "" || echo "_${DSET#depth_}")
OUT=$RL/gwas$SUF
mkdir -p "$OUT"; cd "$OUT"

for f in "$SCRIPTS/make_readlevel_rtab.py" "$SCRIPTS/phenotypes.tsv" "$SCRIPTS/covariates.tsv" \
         "$PY/phylogeny_K.tsv"; do
    [ -s "$f" ] || { echo "missing input: $f"; exit 1; }
done
[ -d "$RL/$DSET" ] || { echo "missing $RL/$DSET - run 08a (or 08c for depth_q0) first"; exit 1; }

echo "===== 1. read-level genotype matrix (depth set: $DSET -> $OUT)"
python3 "$SCRIPTS/make_readlevel_rtab.py" "$RL/$DSET" "$SCRIPTS/covariates.tsv" "$OUT/readlevel"

npat () { sort -u "$1" 2>/dev/null | wc -l; }

echo; echo "===== 2. invariant genes (cannot be association-tested)"
# A gene present in all 257 isolates has zero variance: the LMM is singular and pyseer
# cannot test it. That is not a loss - invariance IS the result for those genes, and it
# is reported here descriptively instead.
awk -F'\t' 'NR>1 {p=0; for(i=2;i<=NF;i++) p+=$i;
    printf "  %-14s present in %d/%d%s\n", $1, p, NF-1,
           (p==NF-1 ? "   <- invariant: present in every isolate, not testable" : "")}' \
    "$OUT/readlevel.Rtab"

echo; echo "===== 3. pyseer LMM on read-level genotypes (variable genes only)"
# 0.01-0.99 drops the invariant genes automatically and keeps everything informative.
run () {   # $1 phenotype column  $2 phenotype file
    local tag="$1_readlevel"
    if ! pyseer --phenotypes "$2" --phenotype-column "$1" --pres "$OUT/readlevel.Rtab" \
           --similarity "$PY/phylogeny_K.tsv" --lmm \
           --min-af 0.01 --max-af 0.99 --cpu 4 \
           --output-patterns "patterns_${tag}.txt" > "pyseer_${tag}.tsv" 2> "pyseer_${tag}.log"; then
        echo "--- $1: pyseer FAILED, its log says:"; sed 's/^/      /' "pyseer_${tag}.log"; echo
        return 0
    fi
    local n thr
    n=$(npat "patterns_${tag}.txt"); thr=$(awk -v n="$n" 'BEGIN{if(n>0) printf "%.3g", 0.05/n; else printf "n/a"}')
    echo "--- $1 (patterns=$n, threshold=$thr)"
    if [ "$(wc -l < "pyseer_${tag}.tsv")" -le 1 ]; then
        echo "      no variants tested; pyseer log:"; sed 's/^/      /' "pyseer_${tag}.log"
    else
        column -t -s$'\t' "pyseer_${tag}.tsv" | cut -c1-150
    fi
    echo
}
run EPTB "$SCRIPTS/phenotypes.tsv"
for c in T3ETH_vs_otherL4222 T3ETH_vs_otherL4 T3ETH_vs_all; do
    [ -s "$SCRIPTS/pheno_$c.tsv" ] && run "$c" "$SCRIPTS/pheno_$c.tsv"
done

echo "===== 4. side by side with the assembly-based results"
printf "%-34s %-14s %-14s\n" variant p_assembly p_readlevel
for v in group_1905 group_1965 group_589 group_1946 group_2137; do
    pa=$(awk -F'\t' -v v="$v" '$1==v{print $4}' "$PY/pyseer_EPTB_genes.tsv" 2>/dev/null | head -1)
    pr=$(awk -F'\t' -v v="$v" '$1==v{print $4}' "pyseer_EPTB_readlevel.tsv" 2>/dev/null | head -1)
    printf "%-34s %-14s %-14s\n" "$v" "${pa:-not tested}" "${pr:-NA}"
done

echo; echo "===== done ($(date)); results: $OUT"