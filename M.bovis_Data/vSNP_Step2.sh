#!/bin/bash
#SBATCH --job-name=vSNP_step2
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Load and activate environment
source ~/.bashrc
conda activate vsnp3

# Define directories
STEP1_OUT="/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples/all_fastqs/vSNP_output/step1_output"
OUT_DIR="$STEP1_OUT/step2_output"
mkdir -p "$OUT_DIR"

# Reference genome used in Step1
REF="/home/ma95362/vsnp3_test_dataset/vsnp_dependencies/Mycobacterium_AF2122/NC_002945v4.fasta"
REF_GBK="/home/ma95362/vsnp3_test_dataset/vsnp_dependencies/Mycobacterium_AF2122/NC_002945v4.gbk"


# Move to Step2 folder
cd "$OUT_DIR"

# Run vSNP Step 2
vsnp3_step2.py -wd "$STEP1_OUT" -o "$OUT_DIR" -t "$REF" -b "$REF_GBK" -a -html_tree

# Deactivate environment
conda deactivate

echo "✅ All vSNP Step2 jobs finished successfully."
