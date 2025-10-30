#!/bin/bash
#SBATCH --job-name=phylo_dist
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# ------------------------------
# 1. Activate Conda environment
# ------------------------------
module load Miniforge3
source activate pyseer-env

# --------------------------
# Define input/output files
# --------------------------
TREEFILE=/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_of_1368/core_tree.treefile         # <-- change this
PYSEER_OUT=/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_of_1368/panaroo/pyseer_output  # <-- change this

# Create output directory if it doesn't exist
mkdir -p ${PYSEER_OUT}

# --------------------------
# Run the phylogenetic distance matrix generation
# --------------------------
echo "Running phylogenetic distance matrix generation..."
phylogeny_distance.py --lmm ${TREEFILE} > ${PYSEER_OUT}/phylogeny_K.tsv

echo "✅ Done! Output saved to: ${PYSEER_OUT}/phylogeny_K.tsv"
# Deactivate conda environment if used
# conda deactivate  # <-- uncomment if you use conda   