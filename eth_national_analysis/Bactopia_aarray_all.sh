#!/bin/bash
#SBATCH --job-name=bactopia_array_full
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=05-00:00:00
#SBATCH --array=1-1398
#SBATCH --output=/scratch/ma95362/eth_national_analysis/logs/%A_%a.out
#SBATCH --error=/scratch/ma95362/eth_national_analysis/logs/%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# -------------------------
# Load module
# -------------------------
module load Bactopia/3.2.0

# -------------------------
# Paths
# -------------------------
SAMPLES=/scratch/ma95362/eth_national_analysis/bactopia_prepare/samples.txt
RESULTS=/scratch/ma95362/eth_national_analysis/bactopia_results
LOGS=/scratch/ma95362/eth_national_analysis/logs

mkdir -p "$RESULTS" "$LOGS"

# -------------------------
# Prepare samples file (skip header & empty lines)
# -------------------------
SAMPLES_CLEAN=$(mktemp)
awk 'NR>1 && NF>0' "$SAMPLES" > "$SAMPLES_CLEAN"  # skip header & empty lines

TOTAL=$(wc -l < "$SAMPLES_CLEAN")
echo "Total samples to process: $TOTAL"

# -------------------------
# Check SLURM_ARRAY_TASK_ID
# -------------------------
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "This script is meant to be run as an array job."
    echo "Please submit with:"
    echo "sbatch --array=1-$TOTAL $0"
    exit 1
fi

# -------------------------
# Get the line for this array task
# -------------------------
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_CLEAN" | tr -d '\r')

SAMPLE=$(echo "$LINE" | cut -f1)
FQ1=$(echo "$LINE" | cut -f5)  # corrected column for R1
FQ2=$(echo "$LINE" | cut -f6)  # corrected column for R2

OUTDIR="$RESULTS/$SAMPLE"
mkdir -p "$OUTDIR"

echo "[$SLURM_ARRAY_TASK_ID] Processing sample: $SAMPLE"
echo "FASTQ R1: $FQ1"
echo "FASTQ R2: $FQ2"
echo "Output directory: $OUTDIR"

# -------------------------
# Skip already-processed samples
# -------------------------
#if [ -f "$OUTDIR/bactopia.done" ]; then
#    echo "Sample $SAMPLE already processed. Skipping."
#    exit 0
#fi
#
#cd "$OUTDIR" || { echo "Failed to cd into $OUTDIR"; exit 1; }
#
# -------------------------
# Run Bactopia
# -------------------------
bactopia \
    --r1 "$FQ1" \
    --r2 "$FQ2" \
    --sample "$SAMPLE" \
    --outdir "$OUTDIR" \
    --species "Mycobacterium tuberculosis" \
    --cpus 16 \
    --genome-size 4400000

# -------------------------
# Mark as done
# -------------------------
touch "$OUTDIR/bactopia.done"
echo "[$SLURM_ARRAY_TASK_ID] Sample $SAMPLE finished."

# -------------------------
# Cleanup
# -------------------------
#rm -f "$SAMPLES_CLEAN"
