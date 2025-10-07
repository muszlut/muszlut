#!/bin/bash
#SBATCH --job-name=Bactopia_Run
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
# Load environment modules
# -----------------------------
module load Bactopia/3.2.0-conda
module load Java/17.0.6

# -----------------------------
# Define directories
# -----------------------------
OUTDIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples"
FOFN="/scratch/ma95362/eth_national_analysis/all_fastq_reads/ForBactopia_FOFN_samples.fofn"

# -----------------------------
# Step 1: Run Bactopia
# -----------------------------
echo "[$(date)] Starting Bactopia run..."
bactopia \
    --samples $FOFN \
    --outdir $OUTDIR \
    --species "Mycobacterium tuberculosis" \
    --genome-size 4410000 \
    --coverage 100 \
    --max_cpus 23 \
    --resume

# -----------------------------
# Step 2: Summarize and generate plots
# -----------------------------
echo "[$(date)] Generating summary reports..."
bactopia summary --bactopia-path $OUTDIR
bactopia plot --bactopia-path $OUTDIR --outdir $OUTDIR/plots

# -----------------------------
# Step 3: Copy key outputs for easy access
# -----------------------------
cp $OUTDIR/bactopia-summary.txt /scratch/ma95362/eth_national_analysis/all_fastq_reads/
cp $OUTDIR/bactopia-exclude.tsv /scratch/ma95362/eth_national_analysis/all_fastq_reads/
cp $OUTDIR/bactopia-report.tsv /scratch/ma95362/eth_national_analysis/all_fastq_reads/

echo "[$(date)] Bactopia run completed successfully."
