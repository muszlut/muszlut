#!/bin/bash
#SBATCH --job-name=mtb_annotate
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/ma95362/Bio_project_publication/logs/log.%j.out
#SBATCH --error=/scratch/ma95362/Bio_project_publication/logs/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu
#
# Collect and annotate all significant pyseer hits (gwas/pyseer + gwas/T3ETH),
# then print the T3-ETH group-defining tables, so one log holds everything to review.
set -euo pipefail

module load Mamba/23.11.0-0
source activate mtb_gwas

PROJ=/scratch/ma95362/Bio_project_publication
SCRIPTS=${SCRIPTS:-$PROJ/gwas/scripts}
GWAS=$PROJ/gwas
PAN=$PROJ/PanPanaroo_merged_results/bactopia-runs/pangenome-20260927-192230/panaroo
GBK=/scratch/ma95362/gbk/ncbi_dataset/data/GCF_000195955.2/genomic.gbk
cd "$GWAS"

for f in "$SCRIPTS/annotate_hits.py" "$PAN/gene_presence_absence.csv" "$GBK" "$SCRIPTS/phenotypes.tsv"; do
    [ -s "$f" ] || { echo "missing input: $f"; exit 1; }
done

echo "===== 1. significant hits, annotated ($(date))"
python3 "$SCRIPTS/annotate_hits.py" "$GWAS" "$PAN" "$GBK" "$SCRIPTS"

echo; echo "===== 2. T3-ETH group-defining variants"
for f in T3ETH/T3ETH_vs_otherL4222.group_defining.tsv T3ETH/T3ETH_vs_otherL4.group_defining.tsv T3ETH/T3ETH_vs_all.group_defining.tsv; do
    echo "--- $f"; [ -s "$f" ] && cat "$f" || echo "(none)"
done

echo; echo "===== 3. T3-ETH vs other L4.2.2.2: 15 most different variants (cols 1-10)"
awk 'NR<=16' T3ETH/T3ETH_vs_otherL4222.all_variants.tsv | cut -f1-10

echo; echo "===== done ($(date)); full table: $GWAS/significant_hits_annotated.tsv"