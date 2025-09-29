#!/bin/bash
#SBATCH --job-name=bactopia_snippy_array
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=23
#SBATCH --mem=240G
#SBATCH --time=05-00:00:00                              # Time limit (HH:MM:SS)
#SBATCH --array=1-1398
#SBATCH --output=/scratch/ma95362/eth_national_analysis/logs/snippy_%A_%a.out
#SBATCH --error=/scratch/ma95362/eth_national_analysis/logs/snippy_%A_%a.err
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
BRESULTS=/scratch/ma95362/eth_national_analysis/bactopia_results
SNIPPY_RESULTS=/scratch/ma95362/eth_national_analysis/snippy_results
REFERENCE=/scratch/ma95362/gbk/ncbi_dataset/data/GCF_000195955.2/genomic.gbk
SNIPPY_LOGS=/scratch/ma95362/eth_national_analysis/logs

mkdir -p "$SNIPPY_RESULTS" "$SNIPPY_LOGS"

# -------------------------
# Prepare samples file (skip header & empty lines)
# -------------------------
SAMPLES_CLEAN=$(mktemp)
awk 'NR>1 && NF>0' "$SAMPLES" > "$SAMPLES_CLEAN"

TOTAL=$(wc -l < "$SAMPLES_CLEAN")
echo "Total samples to process: $TOTAL"

# -------------------------
# Check SLURM_ARRAY_TASK_ID
# -------------------------
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "This script is meant to be run as an array job."
    echo "Please submit with: sbatch --array=1-$TOTAL $0"
    exit 1
fi

# -------------------------
# Get the line for this array task
# -------------------------
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_CLEAN" | tr -d '\r')
SAMPLE=$(echo "$LINE" | cut -f1)

OUTDIR="$SNIPPY_RESULTS/$SAMPLE"
mkdir -p "$OUTDIR"
cd "$OUTDIR" || { echo "Failed to cd into $OUTDIR"; exit 1; }

echo "[$SLURM_ARRAY_TASK_ID] Processing sample: $SAMPLE"

# -------------------------
# Run Snippy via Bactopia
# -------------------------
bactopia tools snippy \
    --bactopia "$BRESULTS" \
    --outdir "$OUTDIR" \
    --sample "$SAMPLE" \
    --cpus 23 \
    --reference "$REFERENCE"

# -------------------------
# Mark as done
# -------------------------
touch "$OUTDIR/snippy.done"
echo "[$SLURM_ARRAY_TASK_ID] Sample $SAMPLE finished."

