#!/bin/bash
#SBATCH --job-name=Snippy_Run
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=23
#SBATCH --mem=240G
#SBATCH --time=06-00:00:00
#SBATCH --output=/scratch/ma95362/eth_national_analysis/logs/snippy_%j.out
#SBATCH --error=/scratch/ma95362/eth_national_analysis/logs/snippy_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Load Bactopia
module load Bactopia/3.2.0

# Directories
OUTDIR="/scratch/ma95362/eth_national_analysis/snippy_results"
REF="/scratch/ma95362/gbk/ncbi_dataset/data/GCF_000195955.2/genomic.gbk"
BRES="/scratch/ma95362/eth_national_analysis/bactopia_results"

# Create output directory
mkdir -p $OUTDIR

# Run Bactopia Snippy workflow on ALL samples
bactopia \
    --wf snippy \
    --reference $REF \
    --bactopia $BRES \
    --outdir $OUTDIR \
    --cpus 23
