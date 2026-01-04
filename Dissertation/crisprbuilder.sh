#!/bin/bash
#SBATCH --job-name=CRISPRbuilder_SRA
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/ma95362/scratch/ncbi_upload.%j.out
#SBATCH --error=/scratch/ma95362/scratch/ncbi_upload.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# -------------------------
# Load environment
# -------------------------
module purge
source ~/.bashrc
conda activate crisprbuilder_tb

# -------------------------
# Paths
# -------------------------
BASE_DIR=/scratch/ma95362/CRISPRbuilder-TB
SEQ_DIR=$BASE_DIR/sequences
SCRIPT=$BASE_DIR/crisprbuilder.py

cd $BASE_DIR || exit 1

echo "Starting CRISPRbuilder-TB (-sra mode)"
echo "Base directory: $BASE_DIR"
echo "Sequences directory: $SEQ_DIR"

# -------------------------
# Run CRISPRbuilder per sample
# -------------------------
for sample_dir in $SEQ_DIR/*; do
    sample=$(basename "$sample_dir")

    fasta=${sample_dir}/${sample}_shuffled.fasta

    if [[ -f "$fasta" ]]; then
        echo "======================================"
        echo "Processing sample: $sample"
        echo "FASTA: $fasta"
        echo "======================================"

        cd "$sample_dir" || exit 1

        python3 $SCRIPT -sra

        cd $BASE_DIR || exit 1
    else
        echo "WARNING: shuffled FASTA not found for $sample"
    fi
done

echo "All samples processed."
