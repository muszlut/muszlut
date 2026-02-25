#!/bin/bash
#SBATCH --job-name=snippy                                      # Job name
#SBATCH --partition=batch                                      # Partition (queue) name
#SBATCH --ntasks=1                                             # Run on a single CPU
#SBATCH --cpus-per-task=32                                     # Number of cores per task
#SBATCH --mem=120gb                                            # Job memory request
#SBATCH --time=07-00:00:00                                     # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err            # Standard error log

#SBATCH --mail-type=END,FAIL                                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                            # Where to send mail	

# -------------------------------
# Step 0: Set Output Directory
# -------------------------------
# Variables
#READS_DIR="/scratch/ma95362/clean_sequences_reads"
OUTDIR="/scratch/ma95362/publication"
REF="/scratch/ma95362/gbk/ncbi_dataset/data/GCF_000195955.2/genomic.gbk"
# Load Bactopia module
module load Bactopia/3.2.0-conda
module load Java/17.0.6
# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"
# Move to working directory
cd $OUTDIR
# Prepare samples for Bactopia
#bactopia prepare \
#    --path "$READS_DIR" \
#    --species "Mycobacterium tuberculosis" \
#    --genome-size 4410000 \
#    > "$OUTDIR/samples.txt"
# Run Bactopia pipeline
#bactopia \
#    --samples "$OUTDIR/samples.txt" \
#    --coverage 100 \
#    --outdir "$OUTDIR" \
#    --max_cpus 16
# Generate summary
#bactopia summary \
#    --bactopia-path "$OUTDIR"
bactopia \
    --wf snippy \
    --reference $REF \
    --bactopia $OUTDIR 
