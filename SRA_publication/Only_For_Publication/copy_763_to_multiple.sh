#!/bin/bash
#SBATCH --job-name=copy_SRR_ERR
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=07-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

set -euo pipefail

# Source and destination
SRC="/scratch/ma95362/Eth_dataset_763/bactopia_run/results"
DEST="/scratch/ma95362/Bio_project_publication/Mtb_Multiple_BioProjects/results_singularity"

# Create destination
mkdir -p "$DEST"

echo "=========================================="
echo "Starting SRR and ERR directory transfer"
echo "Date: $(date)"
echo "=========================================="

# Find SRR and ERR directories
shopt -s nullglob
DIRS=("$SRC"/SRR*/ "$SRC"/ERR*/)

echo "Total directories identified: ${#DIRS[@]}"

# Copy each complete sample directory
for DIR in "${DIRS[@]}"; do

    SAMPLE=$(basename "$DIR")

    echo "Copying sample: $SAMPLE"

    # Remove trailing slash to preserve sample directory name
    rsync -a --info=progress2 "${DIR%/}" "$DEST/"

done

echo "=========================================="
echo "TRANSFER COMPLETED"
echo "Date: $(date)"
echo "=========================================="

# Verification
SOURCE_COUNT=${#DIRS[@]}

DEST_COUNT=$(find "$DEST" -mindepth 1 -maxdepth 1 \
    -type d \( -name "SRR*" -o -name "ERR*" \) | wc -l)

echo "Source sample directories: $SOURCE_COUNT"
echo "Destination sample directories: $DEST_COUNT"

if [ "$SOURCE_COUNT" -eq "$DEST_COUNT" ]; then
    echo "SUCCESS: Directory counts match."
else
    echo "WARNING: Directory counts do not match."
fi

echo "Job completed: $(date)"