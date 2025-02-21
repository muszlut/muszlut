#!/bin/bash
#SBATCH --job-name=Musse_WGS_BacPrep_second                  # Job name
#SBATCH --partition=batch                                    # Partition (queue) name
#SBATCH --ntasks=1                                           # Run on a single CPU
#SBATCH --cpus-per-task=8                                    # Number of cores per task
#SBATCH --mem=40gb                                           # Job memory request
#SBATCH --time=05-00:00:00                                   # Time limit (days-hours:minutes:seconds)
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out        # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err         # Standard error log

#SBATCH --mail-type=END,FAIL                                 # Mail notifications (END, FAIL)
#SBATCH --mail-user=ma95362@uga.edu                          # Email for notifications	

# Define directories
INPUT_DIR="/work/fdqlab/Ethiopia_wgs_mtb_2024/second_run"
OUTDIR="/scratch/ma95362/musse_MGA/fastqs/Second_run"

# Create output directory if it does not exist
if [ ! -d "$OUTDIR" ]; then
    mkdir -p "$OUTDIR"
fi

# Load Bactopia module
module load Bactopia/3.1.0

# Move to the output directory
cd "$OUTDIR"

# Prepare Bactopia samples list with FASTQ patterns
bactopia prepare \
    --path "$INPUT_DIR" \
    --fastq-ext .fq.gz \
    --pe1-pattern .1 \
    --pe2-pattern .2 \
    >$OUTDIR/second_run_samples.txt"

# Run Bactopia with the prepared samples list
bactopia \
    --samples $OUTDIR/second_run_samples.txt" \
    --coverage 100 \
    --outdir $OUTDIR/local_multiple_samples" \
    --max_cpus 4

#bactopia summary \
##   --bactopia-path $OUTDIR/MGA_paired_end_samples

#bactopia search \
#    --query PRJNA1155881
#bactopia \
#    --accessions $OUTDIR/bactopia-accessions.txt \
#    --coverage 100 \
#    --outdir $OUTDIR/ena-multiple-samples \
#    --max_cpus 8