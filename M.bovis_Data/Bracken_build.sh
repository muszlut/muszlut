#!/bin/bash
#SBATCH --job-name=bracken_build          # Job name
#SBATCH --partition=highmem               # High-memory partition
#SBATCH --ntasks=1                        # Single task
#SBATCH --cpus-per-task=8                 # Threads for Bracken
#SBATCH --mem=120gb                        # Memory
#SBATCH --time=12:00:00                   # Time limit (HH:MM:SS)
#SBATCH --output=/scratch/ma95362/bracken_build.%j.out
#SBATCH --error=/scratch/ma95362/bracken_build.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Load Conda and activate Bovisanalyzer environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bovisanalyzer

# Path to Kraken2 DB
KRAKEN2DB=/scratch/ma95362/kraken2_db

# Change directory to Kraken2 DB
cd $KRAKEN2DB || { echo "Failed to cd into $KRAKEN2DB"; exit 1; }

# Build Bracken k-mer distribution
bracken-build \
    -d . \
    -t 8 \
    -k 35 \
    -l 100

echo "Bracken build completed!"
