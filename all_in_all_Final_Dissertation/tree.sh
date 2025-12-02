#!/bin/bash
#SBATCH --job-name=iqtree_core
#SBATCH --partition=batch
#SBATCH --output=iqtree_core_%j.out    # Standard output file
#SBATCH --error=iqtree_core_%j.err     # Standard error file
#SBATCH --ntasks=1                      # Number of tasks (1 for IQ-TREE)
#SBATCH --cpus-per-task=8               # Number of CPU cores for IQ-TREE
#SBATCH --time=12:00:00                 # Max runtime (adjust as needed)
#SBATCH --mem=16G                       # Memory (adjust as needed)
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# ------------------------------
# 1. Activate Conda environment
# ------------------------------
module load Miniforge3
source activate pyseer-env

cd /scratch/ma95362/all_in_all_reads/bactopia_prepare/bactopia-runs/pangenome-20251128-070449
# Run IQ-TREE
iqtree -s core_gene_alignment.aln -pre core_tree -nt 8 -fast -m GTR
