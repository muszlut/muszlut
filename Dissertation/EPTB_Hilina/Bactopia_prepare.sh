#!/bin/bash
#SBATCH --job-name=Bactopia_Mtb
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120gb
#SBATCH --time=03-00:00:00
#SBATCH --output=/scratch/ma95362/log.%j.out
#SBATCH --error=/scratch/ma95362/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Load the Bactopia module
module load Bactopia/3.2.0-conda
module load EDirect/20.5.20231006-GCCcore-12.3.0
# Define directories
OUTDIR="/scratch/ma95362/EPTB_Hilina/Newe"
# Step 1: Create main output directories
mkdir -p $OUTDIR
cd $OUTDIR 
#bactopia search \
#    --query PRJNA1174701 \
#    --force \
#    --genome-size 4400000 \
#    --min-coverage 20
bactopia --accessions bactopia-filtered.txt \
         --outdir $OUTDIR/Bactopia_Run \
         --species "Mycobacterium tuberculosis" \
         --threads 16
