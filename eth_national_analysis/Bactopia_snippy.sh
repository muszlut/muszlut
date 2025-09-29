#!/bin/bash
#SBATCH --job-name=bactopia_snippy_array
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=23
#SBATCH --mem=240G
#SBATCH --time=05-00:00:00
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
SNIPPY_RESULTS=/scratch/ma95362/eth_national_analysis/snippy_results
REFERENCE=/scratch/ma95362/gbk/ncbi_dataset/data/GCF_000195955.2/genomic.gbk

mkdir -p "$SNIPPY_RESULTS"

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
FQ1=$(echo "$LINE" | cut -f5)   # column with R1
FQ2=$(echo "$LINE" | cut -f6)   # column with R2

OUTDIR="$SNIPPY_RESULTS/$SAMPLE"
mkdir -p "$OUTDIR"
cd "$OUTDIR" || { echo "Failed to cd into $OUTDIR"; exit 1; }

echo "[$SLURM_ARRAY_TASK_ID] Processing sample: $SAMPLE"
echo "R1: $FQ1"
echo "R2: $FQ2"
echo "Output directory: $OUTDIR"

# -------------------------
# Run Snippy via main Bactopia pipeline
# -------------------------
bactopia \
    --r1 "$FQ1" \
    --r2 "$FQ2" \
    --sample "$SAMPLE" \
    --outdir "$OUTDIR" \
    --wf snippy \
    --species "Mycobacterium tuberculosis" \
    --reference "$REFERENCE" \
    --cpus 23

# -------------------------
# Mark as done
# -------------------------
touch "$OUTDIR/snippy.done"
echo "[$SLURM_ARRAY_TASK_ID] Sample $SAMPLE finished."


