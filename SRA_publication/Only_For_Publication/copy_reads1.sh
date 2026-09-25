#!/bin/bash
#SBATCH --job-name=copy_E10_E22_E23
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/ma95362/scratch/copy_E10_E22_E23.%j.out
#SBATCH --error=/scratch/ma95362/scratch/copy_E10_E22_E23.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

set -euo pipefail

# Source and destination directories
SRC="/work/fdqlab/Ethiopia_wgs_mtb_2024/first_run"
DEST="/scratch/ma95362/clean_sequences_reads"

# Create destination directory
mkdir -p "$DEST"

# Files to copy
FILES=(
    "E10._R1.fastq.gz"
    "E10._R2.fastq.gz"
    "E22._R1.fastq.gz"
    "E22._R2.fastq.gz"
    "E23._R1.fastq.gz"
    "E23._R2.fastq.gz"
)

# Check that all files exist
for FILE in "${FILES[@]}"; do
    if [ ! -f "$SRC/$FILE" ]; then
        echo "ERROR: File not found: $SRC/$FILE"
        exit 1
    fi
done

# Copy files
for FILE in "${FILES[@]}"; do
    cp -v "$SRC/$FILE" "$DEST/"
done

echo ""
echo "======================================"
echo "Copy completed successfully!"
echo "======================================"

# Verify copied files
ls -lh \
    "$DEST/E10._R1.fastq.gz" \
    "$DEST/E10._R2.fastq.gz" \
    "$DEST/E22._R1.fastq.gz" \
    "$DEST/E22._R2.fastq.gz" \
    "$DEST/E23._R1.fastq.gz" \
    "$DEST/E23._R2.fastq.gz"