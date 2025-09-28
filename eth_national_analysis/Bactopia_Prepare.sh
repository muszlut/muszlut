#!/bin/bash
#SBATCH --job-name=bactopia_prepare                      # Job name
#SBATCH --partition=batch                               # Partition (queue)
#SBATCH --ntasks=1                                      # Single task
#SBATCH --cpus-per-task=16                              # CPUs per task
#SBATCH --mem=120gb                                      # Memory
#SBATCH --time=05-00:00:00                              # Time limit (HH:MM:SS)
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