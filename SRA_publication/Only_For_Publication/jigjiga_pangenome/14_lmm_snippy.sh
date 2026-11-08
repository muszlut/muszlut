#!/bin/bash
#SBATCH --job-name=mtb_lmm_snippy
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
# The like-for-like association test on the step 13 genotypes (coverage at each candidate's
# H37Rv coordinates, taken from snippy's whole-genome alignments).
#
# Step 13's Fisher tests ignore population structure, so they are not comparable with the
# pyseer p-values that produced the original hits. This runs the SAME LMM with the SAME
# kinship matrix and phenotype; only the genotype source differs.
#
# Reference points: the genome-wide threshold for the gene test that selected these variants
# was 0.05/194 patterns = 2.58e-4. A three-gene Bonferroni computed here would be circular,
# because these three were chosen for being top hits.
set -euo pipefail

module load Mamba/23.11.0-0
source activate mtb_gwas

PROJ=/scratch/ma95362/Bio_project_publication
SCRIPTS=${SCRIPTS:-$PROJ/gwas/scripts}
PY=$PROJ/gwas/pyseer
DS=$PROJ/gwas/depth_snippy
OUT=$DS/gwas
mkdir -p "$OUT"; cd "$OUT"

for f in "$SCRIPTS/make_readlevel_rtab.py" "$SCRIPTS/phenotypes.tsv" "$SCRIPTS/covariates.tsv" \
         "$PY/phylogeny_K.tsv"; do
    [ -s "$f" ] || { echo "missing input: $f"; exit 1; }
done
[ -d "$DS/depth" ] || { echo "missing $DS/depth - run 13_depth_from_snippy.sh first"; exit 1; }

echo "===== 1. genotype matrix from the snippy-derived coverage"
python3 "$SCRIPTS/make_readlevel_rtab.py" "$DS/depth" "$SCRIPTS/covariates.tsv" "$OUT/snippy"

echo; echo "===== 2. invariant loci (present everywhere; cannot be association-tested)"
awk -F'\t' 'NR>1 {p=0; for(i=2;i<=NF;i++) p+=$i;
    printf "  %-14s present in %d/%d%s\n", $1, p, NF-1,
           (p==NF-1 ? "   <- invariant, not testable" : "")}' "$OUT/snippy.Rtab"

npat () { sort -u "$1" 2>/dev/null | wc -l; }
echo; echo "===== 3. pyseer LMM (variable loci only)"
run () {   # $1 phenotype column  $2 phenotype file
    local tag="$1_snippy"
    if ! pyseer --phenotypes "$2" --phenotype-column "$1" --pres "$OUT/snippy.Rtab" \
           --similarity "$PY/phylogeny_K.tsv" --lmm \
           --min-af 0.01 --max-af 0.99 --cpu 4 \
           --output-patterns "patterns_${tag}.txt" > "pyseer_${tag}.tsv" 2> "pyseer_${tag}.log"; then
        echo "--- $1: pyseer FAILED:"; sed 's/^/      /' "pyseer_${tag}.log"; echo; return 0
    fi
    echo "--- $1 (patterns=$(npat "patterns_${tag}.txt"))"
    { column -t -s$'\t' "pyseer_${tag}.tsv" 2>/dev/null || cat "pyseer_${tag}.tsv"; } | cut -c1-150
    echo
}
run EPTB "$SCRIPTS/phenotypes.tsv"
for c in T3ETH_vs_otherL4222 T3ETH_vs_otherL4 T3ETH_vs_all; do
    [ -s "$SCRIPTS/pheno_$c.tsv" ] && run "$c" "$SCRIPTS/pheno_$c.tsv"
done

echo "===== 4. EPTB: assembly vs each read-based measurement"
printf "%-16s %-14s %-14s %-14s\n" variant p_assembly p_candonly p_wholegenome
for v in group_1905 group_1946 group_1965 group_589 group_2137; do
    pa=$(awk -F'\t' -v v="$v" '$1==v{print $4}' "$PY/pyseer_EPTB_genes.tsv" 2>/dev/null | head -1)
    pc=$(awk -F'\t' -v v="$v" '$1==v{print $4}' "$PROJ/gwas/readlevel/gwas/pyseer_EPTB_readlevel.tsv" 2>/dev/null | head -1)
    pw=$(awk -F'\t' -v v="$v" '$1==v{print $4}' "pyseer_EPTB_snippy.tsv" 2>/dev/null | head -1)
    printf "%-16s %-14s %-14s %-14s\n" "$v" "${pa:-not tested}" "${pc:-NA}" "${pw:-NA}"
done
echo
echo "  genome-wide threshold for the gene test that selected these: 2.58e-4"

echo; echo "===== done ($(date)); results: $OUT"