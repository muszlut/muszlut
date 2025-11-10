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
module load Bactopia/3.2.0-conda
module load EDirect/20.5.20231006-GCCcore-12.3.0

bactopia search \
    --query PRJNA1174701 \
    --output $OUTDIR/Bactopia_Search \
    --force

# Define directories
OUTDIR="/scratch/ma95362/EPTB_Hilina/Newe"


# Step 1: Create main output directories
mkdir -p $OUTDIR

cd $OUTDIR 
bactopia search \
    --query PRJNA1174701 \
    --output $OUTDIR/Bactopia_Search

#echo "=== Step 2: Preparing samples ==="
#bactopia prepare \
#    --path $OUTDIR/Bactopia_Search \
#    --species "Mycobacterium tuberculosis" \
#    --genome-size 4410000 \
#    --output $OUTDIR/Bactopia_Prepare
# Extract the samples file path
#SAMPLES_FILE="${OUTDIR}/Bactopia_Prepare/samples.txt"

#bactopia \
#    --samples $SAMPLES_FILE \
#    --outdir $RUN_DIR \
#    --coverage 100 \
#    --max_cpus 8
#bactopia summary \
#    --bactopia-path $RUN_DIR

