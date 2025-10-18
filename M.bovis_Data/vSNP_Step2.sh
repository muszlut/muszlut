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

# Directories
STEP1_OUT="/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples/all_fastqs/vSNP_output/step1_output"
OUT_DIR="$STEP1_OUT/step2_output"
mkdir -p "$OUT_DIR"
cd "$STEP1_OUT"
# Reference genome used in Step1 (adjust if needed)
REF="/home/ma95362/vsnp3_test_dataset/vsnp_dependencies/Mycobacterium_AF2122/NC_002945v4.fasta"
REF_GBK="/home/ma95362/vsnp3_test_dataset/vsnp_dependencies/Mycobacterium_AF2122/NC_002945v4.gbk"

# Run vSNP Step2
vsnp3_step2.py \
    -wd "$STEP1_OUT" \
    -o "$OUT_DIR" \
    -t fasta \
    -b "$REF_GBK" \
    -v -dp -html_tree

echo "✅ All vSNP Step2 jobs finished successfully."
