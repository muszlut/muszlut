#!/bin/bash
#SBATCH --job-name=copy_E22_E23
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/ma95362/scratch/copy_E22_E23.%j.out
#SBATCH --error=/scratch/ma95362/scratch/copy_E22_E23.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

set -euo pipefail

# Source and destination directories
SRC="/work/fdqlab/Ethiopia_wgs_mtb_2024/second_run"
DEST="/scratch/ma95362/clean_sequences_reads"

# Create destination directory if needed
mkdir -p "$DEST"

# Check that source directory exists
if [ ! -d "$SRC" ]; then
    echo "ERROR: Source directory does not exist: $SRC"
    exit 1
fi

# Check that all four files exist
for file in E22.1.fq.gz E22.2.fq.gz E23.1.fq.gz E23.2.fq.gz
do
    if [ ! -f "$SRC/$file" ]; then
        echo "ERROR: File not found: $SRC/$file"
        exit 1
    fi
done

# Copy E22 and E23 paired-end reads
cp -v \
    "$SRC/E22.1.fq.gz" \
    "$SRC/E22.2.fq.gz" \
    "$SRC/E23.1.fq.gz" \
    "$SRC/E23.2.fq.gz" \
    "$DEST/"

echo ""
echo "Copy completed successfully."
echo "Destination: $DEST"

# Verify copied files
ls -lh \
    "$DEST/E22.1.fq.gz" \
    "$DEST/E22.2.fq.gz" \
    "$DEST/E23.1.fq.gz" \
    "$DEST/E23.2.fq.gz"