#!/bin/bash
#SBATCH --job-name=mtb_rl_summary
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/ma95362/Bio_project_publication/logs/log.%j.out
#SBATCH --error=/scratch/ma95362/Bio_project_publication/logs/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu
#
# Aggregate the read-level depth tables from 08a and test them against EPTB status.
# Run after the 08a array job has finished.
set -euo pipefail

module load Mamba/23.11.0-0
source activate mtb_gwas

PROJ=/scratch/ma95362/Bio_project_publication
SCRIPTS=${SCRIPTS:-$PROJ/gwas/scripts}
PAN=$PROJ/PanPanaroo_merged_results/bactopia-runs/pangenome-20260927-192230/panaroo
OUT=$PROJ/gwas/readlevel

n=$(ls -1 "$OUT/depth"/*.tsv 2>/dev/null | wc -l)
echo "depth tables: $n"
[ "$n" -gt 0 ] || { echo "no depth tables - run 08a first"; exit 1; }

python3 "$SCRIPTS/readlevel_summary.py" "$OUT/depth" "$SCRIPTS/phenotypes.tsv" "$SCRIPTS/covariates.tsv" \
        "$PAN/gene_presence_absence.Rtab" "$PAN/struct_presence_absence.Rtab" \
        "$OUT/readlevel_summary.tsv" | tee "$OUT/readlevel_report.txt"

echo; echo "===== done ($(date)); table: $OUT/readlevel_summary.tsv"