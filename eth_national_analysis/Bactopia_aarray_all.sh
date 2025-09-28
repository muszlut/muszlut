#!/bin/bash
#SBATCH --job-name=bactopia_array_full
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=05-00:00:00
#SBATCH --array=1-1399
#SBATCH --output=/scratch/ma95362/eth_national_analysis/logs/%A_%a.out
#SBATCH --error=/scratch/ma95362/eth_national_analysis/logs/%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# -------------------------
# Setup
# -------------------------
module load Bactopia/3.2.0

# Paths
SAMPLES=/scratch/ma95362/eth_national_analysis/bactopia_prepare/samples.txt
RESULTS=/scratch/ma95362/eth_national_analysis/bactopia_results
LOGS=/scratch/ma95362/eth_national_analysis/logs

mkdir -p "$RESULTS" "$LOGS"

# -------------------------
# Get the line for this array task
# Format of samples.txt: SAMPLE<TAB>R1<TAB>R2
# -------------------------
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES" | tr -d '\r')

SAMPLE=$(echo "$LINE" | cut -f1)
FQ1=$(echo "$LINE" | cut -f2)
FQ2=$(echo "$LINE" | cut -f3)

OUTDIR="$RESULTS/$SAMPLE"
mkdir -p "$OUTDIR"

# Debugging info
echo "[$SLURM_ARRAY_TASK_ID] Processing sample: $SAMPLE"
echo "FASTQ R1: $FQ1"
echo "FASTQ R2: $FQ2"
echo "Output directory: $OUTDIR"

# -------------------------
# Skip already-processed samples
# -------------------------
if [ -f "$OUTDIR/bactopia.done" ]; then
    echo "Sample $SAMPLE already processed. Skipping."
    exit 0
fi

# -------------------------
# Run Bactopia
# -------------------------
bactopia \
    --R1 "$FQ1" \
    --R2 "$FQ2" \
    --sample "$SAMPLE" \
    --outdir "$OUTDIR" \
    --species "Mycobacterium tuberculosis" \
    --genome-size 4400000 \
    --cpus 8

# -------------------------
# Mark as done
# -------------------------
touch "$OUTDIR/bactopia.done"
echo "[$SLURM_ARRAY_TASK_ID] Sample $SAMPLE finished."
