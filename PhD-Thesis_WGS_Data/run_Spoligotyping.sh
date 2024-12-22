#!/bin/bash
#SBATCH --job-name=Musse_Spoligotyping                     # Job name
#SBATCH --partition=batch                                  # Partition (queue) name
#SBATCH --ntasks=1                                         # Run on a single CPU
#SBATCH --cpus-per-task=8                                  # Number of cores per task
#SBATCH --mem=16gb                                         # Job memory request
#SBATCH --time=03-00:00:00                                 # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out       # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err        # Standard error log
#SBATCH --mail-type=END,FAIL                               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                        # Where to send mail 

# Variables
OUTDIR="/scratch/ma95362/musse_MGA/fastqs"
FASTQ_DIR="/scratch/ma95362/musse_MGA/fastqs/MGA_paired_end_samples"               # Directory containing FASTQ files

# Create output directory if it doesn't exist
if [ ! -d $OUTDIR ]; then
    mkdir -p $OUTDIR
fi

# Load the Bactopia module
module load Bactopia/3.1.0

#activate conda envirnoment (should be version 6.3.0)

source activate spotyping-env

# Run the Bactopia spoligotyping workflow
bactopia \
    --wf spotyping \
    --exclude $OUTDIR/bactopia-exclude.tsv \
    --bactopia $FASTQ_DIR \
    --outdir $OUTDIR
