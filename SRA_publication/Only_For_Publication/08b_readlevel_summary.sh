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

# Which depth set to summarise:
#   (no argument) -> depth      MAPQ >= 20, written by 08a and 08c
#   depth_q0      -> depth_q0   unfiltered, written by 08c only
DSET=${1:-depth}
SUF=$([ "$DSET" = "depth" ] && echo "" || echo "_${DSET#depth_}")

n=$(ls -1 "$OUT/$DSET"/*.tsv 2>/dev/null | wc -l)
echo "depth set: $DSET   tables: $n"
[ "$n" -gt 0 ] || { echo "no depth tables in $OUT/$DSET - run 08a (or 08c for depth_q0) first"; exit 1; }

python3 "$SCRIPTS/readlevel_summary.py" "$OUT/$DSET" "$SCRIPTS/phenotypes.tsv" "$SCRIPTS/covariates.tsv" \
        "$PAN/gene_presence_absence.Rtab" "$PAN/struct_presence_absence.Rtab" \
        "$OUT/readlevel_summary${SUF}.tsv" | tee "$OUT/readlevel_report${SUF}.txt"

echo; echo "===== done ($(date)); table: $OUT/readlevel_summary${SUF}.tsv"