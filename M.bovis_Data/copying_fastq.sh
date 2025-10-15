#!/bin/bash
#SBATCH --job-name=copy_fastqs
#SBATCH --output=copy_fastqs_%j.log
#SBATCH --error=copy_fastqs_%j.err
#SBATCH --time=02:00:00          # adjust based on total data size
#SBATCH --mem=8G                 # memory for the job
#SBATCH --cpus-per-task=2        # optional
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Load any necessary modules (none needed for basic cp)
# module load ...

# Set target directory
TARGET_DIR="/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples/all_fastqs"
mkdir -p $TARGET_DIR
cd $TARGET_DIR
# Path to CSV
SAMPLESHEET="/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples/samplesheet.csv"

# Loop through CSV and copy files
tail -n +2 "$SAMPLESHEET" | while IFS=',' read sample fastq1 fastq2; do
    cp "$fastq1" "$TARGET_DIR/"
    cp "$fastq2" "$TARGET_DIR/"
done

