#!/bin/bash
#SBATCH --job-name=copy_SRR_dirs
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=07-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

set -euo pipefail

# Define source and destination
SRC="/scratch/ma95362/Bio_project_publication"
DEST="/scratch/ma95362/Bio_project_publication/Mtb_Multiple_BioProjects/results_singularity"

# Create destination directory
mkdir -p "$DEST"

echo "Copy started: $(date)"

# Copy all SRR directories using rsync
for DIR in "$SRC"/SRR*/; do
    if [ -d "$DIR" ]; then
        echo "Copying $(basename "$DIR") ..."
        rsync -a "$DIR" "$DEST/"

    fi
done

echo "Copy completed: $(date)"

# Verify directory counts
echo "Source SRR directories:"
find "$SRC" -mindepth 1 -maxdepth 1 -type d -name "SRR*" | wc -l

echo "Destination SRR directories:"
find "$DEST" -mindepth 1 -maxdepth 1 -type d -name "SRR*" | wc -l

echo "Job finished."
