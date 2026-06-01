#!/bin/bash
#SBATCH --job-name=Bactopia_prepare
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120gb
#SBATCH --time=07-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu
# -------------------------------
# Variables
# -------------------------------
READS_DIR="/scratch/ma95362/clean_sequences_reads"
OUTDIR="/scratch/ma95362/SRA_publication"

# -------------------------------
# Load Bactopia
# -------------------------------
module load Bactopia/4.0.0-conda

# -------------------------------
# Use lscratch for Bactopia cache
# -------------------------------
#export BACTOPIA_CACHE=/lscratch/ma95362/bactopia-cache
#mkdir -p $BACTOPIA_CACHE/{conda,datasets,singularity}

# -------------------------------
# Create output directory
# -------------------------------
mkdir -p "$OUTDIR"

# -------------------------------
# Move to output directory
# -------------------------------
cd "$OUTDIR"

# -------------------------------
# Prepare samples
# -------------------------------
bactopia prepare \
    --path "$READS_DIR" \
    --species "Mycobacterium tuberculosis" \
    --genome-size 4410000 \
    > "$OUTDIR/samples.txt"

# -------------------------------
# Run Bactopia
# -------------------------------
bactopia \
    --samples "$OUTDIR/samples.txt" \
    --coverage 100 \
    --outdir "$OUTDIR" \
    --max_cpus 16
# Generate summary
#bactopia summary \
#    --bactopia-path "$OUTDIR"
#bactopia \
#    --wf snippy \
#    --reference $REF \
#    --bactopia $OUTDIR 
#bactopia \
#    --wf pangenome \
#    --bactopia $OUTDIR