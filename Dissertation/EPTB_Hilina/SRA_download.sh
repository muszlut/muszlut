#!/bin/bash
#SBATCH --job-name=download_PRJNA1104194
#SBATCH --partition=batch
#SBATCH --ntasks=1                                       # Single task
#SBATCH --cpus-per-task=16                               # CPUs per task
#SBATCH --mem=120GB                                      # Memory
#SBATCH --time=05-00:00:00                               # Time limit (HH:MM:SS)
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out    # STDOUT log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err     # STDERR log
#SBATCH --mail-type=END,FAIL                             # Mail events
#SBATCH --mail-user=ma95362@uga.edu                      # Your email

module load sratoolkit/3.0.0
module load parallel

# Directories
DEST=/scratch/ma95362/EPTB_Hilina/
mkdir -p $DEST
cd $DEST 

# Step 1: Download all SRR files
prefetch --option-file SRR_Acc_List.txt

# Step 2: Convert SRA to FASTQ (gzip-compressed) in parallel
cat SRR_Acc_List.txt | parallel -j 8 "fasterq-dump --split-files --gzip {} -O $DEST"

echo "All sequences downloaded and converted to FASTQ in $DEST"
