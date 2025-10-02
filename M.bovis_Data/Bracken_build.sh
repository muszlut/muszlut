#!/bin/bash
#SBATCH --job-name=bracken_build
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=02-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Load Conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bovisanalyzer

# Path
KRAKEN2DB=/scratch/ma95362/kraken2_db

# Bracken k-mer distribution build
bracken-build \
    -d $KRAKEN2DB \
    -t 16 \
    -k 35 \
    -l 100

echo "Bracken build completed!"