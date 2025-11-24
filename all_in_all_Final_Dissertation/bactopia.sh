#!/bin/bash
#SBATCH --job-name=All_Pangenome
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/ma95362/all_in_all_reads/bactopia_prepare/logs/log.%j.out
#SBATCH --error=/scratch/ma95362/all_in_all_reads/bactopia_prepare/logs/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# -------------------------------
# Step 0: Set Output Directory
# -------------------------------
# Variables
READS_DIR="/scratch/ma95362/all_in_all_reads"
OUTDIR="/scratch/ma95362/all_in_all_reads/bactopia_prepare"
# Load Bactopia module
module load Bactopia/3.2.0-conda
# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"
# Move to working directory
cd $OUTDIR
# Prepare samples for Bactopia
#bactopia prepare \
#    --path "$READS_DIR" \
#    --species "Mycobacterium tuberculosis" \
#    --genome-size 4410000 \
#    > "$OUTDIR/samples.txt"
# Run Bactopia pipeline
bactopia \
    --samples "$OUTDIR/samples.txt" \
    --coverage 100 \
    --outdir "$OUTDIR" \
    --max_cpus 8
# Generate summary
bactopia summary \
    --bactopia-path "$OUTDIR"
#bactopia \
#    --wf tbprofiler \
#    --bactopia $OUTDIR \
#    --max_cpus 32