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
# Set output directory
# -----------------------------
OUTDIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads"
mkdir -p $OUTDIR
cd $OUTDIR

# -----------------------------
# Prepare sample list (only if not exists)
# -----------------------------
if [ ! -f $OUTDIR/ETH_samples.txt ]; then
    bactopia prepare \
        --path $OUTDIR \
        --species "Mycobacterium tuberculosis" \
        --genome-size 4410000 \
        > $OUTDIR/ETH_samples.txt
fi

# -----------------------------
# Run Bactopia on samples with resume and retries
# -----------------------------
bactopia \
    --samples $OUTDIR/ETH_samples.txt \
    --coverage 100 \
    --outdir $OUTDIR/ETH_paired_end_samples \
    --max_cpus 23 \
    --resume

# -----------------------------
# Generate summary and plots
# -----------------------------
bactopia summary \
    --bactopia-path $OUTDIR/ETH_paired_end_samples

bactopia plot \
    --bactopia-path $OUTDIR/ETH_paired_end_samples \
    --outdir $OUTDIR/ETH_paired_end_samples/plots
