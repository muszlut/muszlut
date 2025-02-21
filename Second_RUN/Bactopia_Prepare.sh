#!/bin/bash
#SBATCH --job-name=Musse_WGS_BacPrep_second_run      # Job name for second run
#SBATCH --partition=batch                             # Partition (queue) name
#SBATCH --ntasks=1                                    # Run on a single CPU
#SBATCH --cpus-per-task=8                             # Number of cores per task
#SBATCH --mem=40gb                                    # Job memory request
#SBATCH --time=05-00:00:00                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log_second_run.%j.out  # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log_second_run.%j.err   # Standard error log
#SBATCH --mail-type=END,FAIL                         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                  # Where to send mail

# Set output directory variable for the second run
OUTDIR="/scratch/ma95362/musse_MGA/fastqs_second_run"

# Create the output directory if it doesn't exist
if [ ! -d "$OUTDIR" ]; then
    mkdir -p "$OUTDIR"
fi

# Load the Bactopia module
module load Bactopia/3.1.0

# Change to the output directory
cd "$OUTDIR"

# Prepare sample list from the second run files, adjusting file extension and separator
bactopia prepare \
    --path /work/fdqlab/Ethiopia_wgs_mtb_2024/second_run \
    --species "Mycobacterium tuberculosis" \
    --genome-size 4410000 \
    --fastq-ext .fq.gz \
    --fastq-separator . \
    > "$OUTDIR/MGA_samples_second_run.txt"

# Run Bactopia analysis on the prepared sample list
bactopia \
    --samples "$OUTDIR/MGA_samples_second_run.txt" \
    --coverage 100 \
    --outdir "$OUTDIR/MGA_paired_end_samples_second_run" \
    --max_cpus 4

# Generate a summary of the Bactopia run
bactopia summary \
    --bactopia-path "$OUTDIR/MGA_paired_end_samples_second_run"
