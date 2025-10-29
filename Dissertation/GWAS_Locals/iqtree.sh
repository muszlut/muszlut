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
# 1. module load IQ-TREE
# ------------------------------
module load IQ-TREE/2.3.6-gompi-2023a

# --------------------------
# Define input/output files
# --------------------------
ALIGNMENT=/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_of_1368/core_gene_alignment_filtered.aln
OUTPUT_DIR=/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_of_1368
# --------------------------
# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}
# --------------------------
# Run IQ-TREE for phylogenetic tree construction
# --------------------------
echo "Running IQ-TREE on alignment: $(basename $ALIGNMENT)..."
iqtree2 -s ${ALIGNMENT} -pre core_tree -nt 8 -fast -m GTR

echo "✅ IQ-TREE analysis completed successfully."
echo "Output files saved to: ${OUTPUT_DIR}"
# ------------------------------
# 4. End of script
# ------------------------------