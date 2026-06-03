#!/bin/bash
#SBATCH --job-name=Bactopia_search
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120gb
#SBATCH --time=07-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu
# -------------------------------
# Variables
# -------------------------------
OUTDIR="/scratch/ma95362/SRA_publication/ena_PRJNA1392122"

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
#bactopia search \
#    --query PRJNA1392122
bactopia \
    --accessions bactopia-accessions.txt \
    --coverage 100 \
    --outdir ena-multiple-samples \
    --max_cpus 16