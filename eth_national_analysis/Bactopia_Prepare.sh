#!/bin/bash
#SBATCH --job-name=bactopia_prepare
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out    # STDOUT log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err     # STDERR log
#SBATCH --mail-type=END,FAIL                             # Mail events
#SBATCH --mail-user=ma95362@uga.edu                      # Your email


module load Bactopia/3.2.0
# Define paths
OUTDIR=/scratch/ma95362/eth_national_analysis/bactopia_prepare

# Make output directory if it doesn't exist
mkdir -p $OUTDIR
cd $OUTDIR


bactopia prepare \
    --path /scratch/ma95362/eth_national_analysis/all_fastq_reads \
    --species "Mycobacterium tuberculosis" \
    --genome-size 4400000 \
    > /scratch/ma95362/eth_national_analysis/bactopia_prepare/samples.txt