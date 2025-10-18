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
# Load excluded sample names into memory
EXCLUDE_FILE="/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples/all_fastqs/exclude_list.txt"
Step1="/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples/all_fastqs/vSNP_output"/step1_output
mkdir -p "$OUT_DIR"
cd "$STEP1_OUT"
# Reference genome used in Step1 (adjust if needed)
REF="/home/ma95362/vsnp3_test_dataset/vsnp_dependencies/Mycobacterium_AF2122/NC_002945v4.fasta"
REF_GBK="/home/ma95362/vsnp3_test_dataset/vsnp_dependencies/Mycobacterium_AF2122/NC_002945v4.gbk"

# Run vSNP Step2
#Find all vcf files in the step1 output directory and move them in step2 directory
find $Step1 -type f -name "*_zc.vcf" -exec cp {} $Step2 \;
#move to step2 folder
cd $OUT_DIR
# Run vSNP step 2 to combine SNPs and build tree
vsnp3_step2.py -a -t $REF -remove_by_name $EXCLUDE_FILE
