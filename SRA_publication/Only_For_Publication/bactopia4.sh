#!/bin/bash
#SBATCH --job-name=bactopia_test
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Load Mamba
module load Mamba/23.11.0-0

# Activate Bactopia environment
source activate bactopia4

# Set output directory
OUTDIR="/scratch/ma95362/Bactopia4_test"

# Create output directory if it does not exist
mkdir -p "$OUTDIR"

# Move to output directory
cd "$OUTDIR" || exit 1

# Run Bactopia test
bactopia -profile test,standard