#!/bin/bash
#SBATCH --job-name=copy_reads_Eth_763
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --time=07-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Exit if any command fails
set -e

# Define source and destination directories
SOURCE="/work/fdqlab/Ethiopia_dataset_763"
DEST="/scratch/ma95362/Eth_dataset_763"

# Create destination directory
mkdir -p "$DEST"

echo "Starting file transfer..."
echo "Source: $SOURCE"
echo "Destination: $DEST"
date

# Copy files and directories while preserving attributes
rsync -avh --progress "$SOURCE/" "$DEST/"

echo "File transfer completed successfully!"
date

# Check transferred size
du -sh "$SOURCE"
du -sh "$DEST"