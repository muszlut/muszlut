#!/bin/bash
#SBATCH --job-name=copy_reads_parallel
#SBATCH --partition=batch                                # Partition (queue)
#SBATCH --ntasks=1                                       # Single task
#SBATCH --cpus-per-task=16                               # CPUs per task
#SBATCH --mem=40gb                                       # Memory
#SBATCH --time=05-00:00:00                               # Time limit (HH:MM:SS)
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out    # STDOUT log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err     # STDERR log
#SBATCH --mail-type=END,FAIL                             # Mail events
#SBATCH --mail-user=ma95362@uga.edu                      # Your email

# Define paths
SRC_DIR=/scratch/ma95362/eth_national_analysis/all_fastq_reads
DEST_DIR=/scratch/ma95362/ggcaller_reads

# Create destination if not exist
mkdir -p $DEST_DIR

# Copy all fastq.gz files
cp ${SRC_DIR}/*.fastq.gz ${DEST_DIR}/
