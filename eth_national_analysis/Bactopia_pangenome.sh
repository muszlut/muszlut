#!/bin/bash
#SBATCH --job-name=bactopia_pangenome_array
#SBATCH --partition=highmem_p           # high memory partition for clustering
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32             # parallelize clustering
#SBATCH --mem=240G                     # enough memory for 1398 genomes
#SBATCH --time=7-00:00:00
#SBATCH --array=1-1398                 # one job per sample
#SBATCH --output=/scratch/ma95362/eth_national_analysis/logs/pangenome_%A_%a.out
#SBATCH --error=/scratch/ma95362/eth_national_analysis/logs/pangenome_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# -------------------------
# Load Bactopia
# -------------------------
module load Bactopia/3.2.0

# -------------------------
# Paths
# -------------------------
SAMPLES=/scratch/ma95362/eth_national_analysis/bactopia_prepare/samples.txt
BRESULTS=/scratch/ma95362/eth_national_analysis/bactopia_results
PANRESULTS=/scratch/ma95362/eth_national_analysis/pangenome_results
PAN_LOGS=/scratch/ma95362/eth_national_analysis/logs

mkdir -p "$PANRESULTS" "$PAN_LOGS"

# -------------------------
# Prepare samples file
# -------------------------
SAMPLES_CLEAN=$(mktemp)
awk 'NR>1 && NF>0' "$SAMPLES" > "$SAMPLES_CLEAN"
TOTAL=$(wc -l < "$SAMPLES_CLEAN")
echo "Total samples: $TOTAL"

# -------------------------
# Check SLURM_ARRAY_TASK_ID
# -------------------------
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "This script is meant to be run as an array job."
    echo "Submit with: sbatch --array=1-$TOTAL $0"
    exit 1
fi

# -------------------------
# Per-sample processing (gene extraction, QC)
# -------------------------
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_CLEAN" | tr -d '\r')
SAMPLE=$(echo "$LINE" | cut -f1)
OUTDIR="$PANRESULTS/$SAMPLE"
mkdir -p "$OUTDIR"
cd "$OUTDIR" || { echo "Failed to cd into $OUTDIR"; exit 1; }

echo "[$SLURM_ARRAY_TASK_ID] Processing sample: $SAMPLE"

bactopia tools pangenome \
    --bactopia "$BRESULTS" \
    --sample "$SAMPLE" \
    --outdir "$OUTDIR" \
    --cpus 16

touch "$OUTDIR/pangenome.done"
echo "[$SLURM_ARRAY_TASK_ID] Sample $SAMPLE done."

# -------------------------
# Optional: Final full pangenome clustering
# Only run this once after all samples are processed
# -------------------------
if [ "$SLURM_ARRAY_TASK_ID" -eq "$TOTAL" ]; then
    echo "All samples done. Running final clustering step..."
    FINAL_OUT="$PANRESULTS/final_pangenome"
    mkdir -p "$FINAL_OUT"
    bactopia tools pangenome \
        --bactopia "$BRESULTS" \
        --outdir "$FINAL_OUT" \
        --cpus 32
    echo "Final pangenome clustering done."
fi
