#!/bin/bash
#SBATCH --job-name=mtb_gwas_all
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/ma95362/Bio_project_publication/logs/log.%j.out
#SBATCH --error=/scratch/ma95362/Bio_project_publication/logs/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu
#
# Launcher kept in the home directory; all step scripts and data files are in
#   /scratch/ma95362/Bio_project_publication/gwas/scripts
# Usage:
#   sbatch run-all_scripts.sh               (tree already built)
#   sbatch run-all_scripts.sh --with-tree   (rebuild the masked tree first)
set -euo pipefail

PROJ=/scratch/ma95362/Bio_project_publication
export SCRIPTS=$PROJ/gwas/scripts          # step scripts read their inputs from here
cd "$SCRIPTS"

if [ "${1:-}" = "--with-tree" ]; then
    echo "===== 02_masked_tree.sh  ($(date))"
    bash "$SCRIPTS/02_masked_tree.sh"
fi

echo "===== 03_pyseer_gwas.sh  ($(date))"
bash "$SCRIPTS/03_pyseer_gwas.sh"

echo "===== 04_T3ETH_compare.sh ($(date))"
bash "$SCRIPTS/04_T3ETH_compare.sh"

echo "===== all done ($(date))"