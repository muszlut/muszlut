#!/bin/bash
#SBATCH --job-name=Bactopia_Mtb
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120gb
#SBATCH --time=03-00:00:00
#SBATCH --output=/scratch/ma95362/log.%j.out
#SBATCH --error=/scratch/ma95362/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Load the Bactopia module
module load Bactopia/3.1.0

# Define directories
OUTDIR="/scratch/ma95362/EPTB_Hilina/Newe"
SEARCH_DIR="${OUTDIR}/Bactopia_Search"
PREPARE_DIR="${OUTDIR}/Bactopia_Prepare"
RUN_DIR="${OUTDIR}/Bactopia_Run"

# Step 1: Create main output directories
mkdir -p $SEARCH_DIR $PREPARE_DIR $RUN_DIR

echo "=== Step 1: Running bactopia search ==="
bactopia search \
    --query PRJNA1174701 \
    --output $SEARCH_DIR

echo "=== Step 2: Preparing samples ==="
bactopia prepare \
    --path $SEARCH_DIR \
    --species "Mycobacterium tuberculosis" \
    --genome-size 4410000 \
    --output $PREPARE_DIR

# Extract the samples file path
SAMPLES_FILE="${PREPARE_DIR}/samples.txt"

echo "=== Step 3: Running bactopia main analysis ==="
bactopia \
    --samples $SAMPLES_FILE \
    --outdir $RUN_DIR \
    --coverage 100 \
    --max_cpus 8

echo "=== Step 4: Generating summary report ==="
bactopia summary \
    --bactopia-path $RUN_DIR

echo "=== All steps completed successfully! ==="
