#!/bin/bash
#SBATCH --job-name=DBactopia_Run
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=23
#SBATCH --mem=240G
#SBATCH --time=07-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# -----------------------------
# Load required modules
# -----------------------------
module load Bactopia/3.2.0-conda
module load Java/17.0.6

# -----------------------------
# Set directories and FOFN
# -----------------------------
OUTDIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads"
FOFN="$OUTDIR/TBprofiler_FOFN_samples.fofn"
RESULTS_DIR="$OUTDIR/ETH_paired_end_samples"

mkdir -p $RESULTS_DIR
cd $OUTDIR

# -----------------------------
# Clean old reports to avoid conflicts
# -----------------------------
OLD_REPORTS=(
    "$RESULTS_DIR/bactopia-report.tsv"
    "$RESULTS_DIR/bactopia-exclude.tsv"
    "$RESULTS_DIR/bactopia-summary.txt"
)

for file in "${OLD_REPORTS[@]}"; do
    if [ -f "$file" ]; then
        echo "Removing old report: $file"
        rm "$file"
    fi
done

# -----------------------------
# Run Bactopia using the FOFN
# -----------------------------
bactopia \
    --samples $FOFN \
    --coverage 100 \
    --outdir $RESULTS_DIR \
    --max_cpus 23 \
    --resume

# -----------------------------
# Generate summary and plots
# -----------------------------
bactopia summary \
    --bactopia-path $RESULTS_DIR

bactopia plot \
    --bactopia-path $RESULTS_DIR \
    --outdir $RESULTS_DIR/plots
