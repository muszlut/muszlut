#!/bin/bash
#SBATCH --job-name=copy_SRR2882_files
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/ma95362/EPTB_Hilina/new_project_logs/Bactopia_prepare/log.%j.out
#SBATCH --error=/scratch/ma95362/EPTB_Hilina/new_project_logs/Bactopia_prepare/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

SRC_DIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads"
DEST_DIR="/scratch/ma95362/all_in_all_reads"

echo "Copying all SRR2882* fastq.gz files..."
mkdir -p "$DEST_DIR"

cp ${SRC_DIR}/SRR2882.fastq.gz "$DEST_DIR"/

echo "Copy completed."

