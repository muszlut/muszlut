#!/bin/bash
#SBATCH --job-name=iqtree_core
#SBATCH --partition=batch
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --mem=160G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# ------------------------------
# 1. Activate Conda environment
# ------------------------------
source $(conda info --base)/etc/profile.d/conda.sh
conda activate panaroo-env

# ------------------------------
# 2. Change to working directory
# ------------------------------
cd /scratch/ma95362/all_in_all_reads/bactopia_prepare/bactopia-runs/pangenome-20251128-070449

# ------------------------------
# 3. Run IQ-TREE
# ------------------------------
iqtree -s core_gene_alignment_filtered.aln -pre core_tree -nt 8 -fast -m GTR -safe
