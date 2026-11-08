#!/bin/bash
#SBATCH --job-name=mtb_eptb_cov_tree
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/ma95362/Bio_project_publication/logs/log.%j.out
#SBATCH --error=/scratch/ma95362/Bio_project_publication/logs/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu
#
# (a) EPTB vs PTB rerun with sequencing coverage and contig number as covariates
#     -> tests whether the EPTB gene/structural hits are batch (coverage) artefacts
# (b) phylogeny figure with T3-ETH marker tracks -> gwas/figures/
set -euo pipefail

module load Mamba/23.11.0-0
source activate mtb_gwas

PROJ=/scratch/ma95362/Bio_project_publication
SCRIPTS=${SCRIPTS:-$PROJ/gwas/scripts}
PAN=$PROJ/PanPanaroo_merged_results/bactopia-runs/pangenome-20260927-192230/panaroo
PY=$PROJ/gwas/pyseer                       # results of step 03 (kinship, masked VCF, uncovaried runs)
TREE=$PROJ/gwas/tree/mtb257_masked.treefile.midpoint_rooted.nwk
OUT=$PROJ/gwas/eptb_covariates
FIG=$PROJ/gwas/figures
mkdir -p "$OUT" "$FIG"

for f in "$SCRIPTS/phenotypes.tsv" "$SCRIPTS/covariates.tsv" "$SCRIPTS/isolate_annotations.tsv" "$SCRIPTS/plot_tree.py" \
         "$PY/phylogeny_K.tsv" "$PY/core_snps.masked.vcf.gz" "$TREE" "$PROJ/gwas/T3ETH/snps_matrix.tsv" \
         "$PAN/gene_presence_absence.Rtab" "$PAN/struct_presence_absence.Rtab"; do
    [ -s "$f" ] || { echo "missing input: $f"; exit 1; }
done

npat () { sort -u "$1" | wc -l; }
thresh () { awk -v n="$(npat "$1")" 'BEGIN{printf "%.3g", 0.05/n}'; }

echo "===== (a) EPTB with covariates: coverage + contigs ($(date))"
cd "$OUT"; : > summary.txt
run () {   # $1 label  $2 pyseer input flag  $3 file
    local tag="EPTB_$1_cov"
    pyseer --phenotypes "$SCRIPTS/phenotypes.tsv" --phenotype-column EPTB $2 "$3" \
           --similarity "$PY/phylogeny_K.tsv" --lmm \
           --covariates "$SCRIPTS/covariates.tsv" --use-covariates 2q 3q \
           --min-af 0.05 --max-af 0.95 --cpu 8 \
           --output-patterns "patterns_${tag}.txt" > "pyseer_${tag}.tsv" 2> "pyseer_${tag}.log"
    local thr; thr=$(thresh "patterns_${tag}.txt")
    local nsig; nsig=$(awk -F'\t' -v t="$thr" 'NR>1 && $4!="" && $4<t' "pyseer_${tag}.tsv" | wc -l)
    printf "%-22s patterns=%-6s threshold=%-10s significant=%s\n" "$tag" \
           "$(npat "patterns_${tag}.txt")" "$thr" "$nsig" | tee -a summary.txt
    sort -t$'\t' -g -k4,4 <(tail -n +2 "pyseer_${tag}.tsv") > "sorted_${tag}.tsv"
    # side-by-side: top 10 without covariates and their p-value with covariates
    local base="$PY/pyseer_EPTB_$1.tsv"
    if [ -s "$base" ]; then
        echo "-- $1: top 10 without covariates vs with covariates (lrt-pvalue, beta)"
        printf "%-45s %-12s %-10s %-12s %-10s\n" variant p_nocov beta_nocov p_cov beta_cov
        sort -t$'\t' -g -k4,4 <(tail -n +2 "$base") | awk 'NR<=10' > tmp_top.tsv
        awk -F'\t' 'NR==FNR{p[$1]=$4; b[$1]=$5; next} {printf "%-45s %-12s %-10s %-12s %-10s\n", substr($1,1,45), $4, $5, ($1 in p?p[$1]:"NA"), ($1 in b?b[$1]:"NA")}' \
            "pyseer_${tag}.tsv" tmp_top.tsv
        rm -f tmp_top.tsv
    fi
    echo "-- $1: top 5 with covariates"; awk 'NR<=5' "sorted_${tag}.tsv" | cut -f1-5; echo
}
run genes  --pres "$PAN/gene_presence_absence.Rtab"
run struct --pres "$PAN/struct_presence_absence.Rtab"
run snps   --vcf  "$PY/core_snps.masked.vcf.gz"

echo "===== (b) phylogeny figure ($(date))"
python3 "$SCRIPTS/plot_tree.py" "$TREE" "$SCRIPTS/isolate_annotations.tsv" "$PROJ/gwas/T3ETH/snps_matrix.tsv" \
        "$PAN/struct_presence_absence.Rtab" "$FIG/T3ETH_phylogeny"

echo; echo "===== done ($(date))"
echo "covariate results: $OUT   figure: $FIG/T3ETH_phylogeny.png / .pdf"