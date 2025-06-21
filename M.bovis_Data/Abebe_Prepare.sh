#!/bin/bash
#SBATCH --job-name=Musse_Mbovis_Bactopia                     # Job name
#SBATCH --partition=batch                                    # Partition (queue) name
#SBATCH --ntasks=1                                           # Run on a single CPU
#SBATCH --cpus-per-task=8                                    # Number of cores per task
#SBATCH --mem=40gb                                           # Job memory request
#SBATCH --time=05-00:00:00                                   # Time limit d-hh:mm:ss
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out         # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err          # Standard error log
#SBATCH --mail-type=END,FAIL                                 # Mail notifications
#SBATCH --mail-user=ma95362@uga.edu                          # Your email address

# Variables
READS_DIR="/scratch/ma95362/ETH_bovis_Sequence/reads/64_SH_Sequence_data/raw/sub_raw"
OUTDIR="/scratch/ma95362/ETH_bovis_Sequence/Abebe_samples"

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Load Bactopia module
module load Bactopia/3.1.0

# Move to working directory
cd $OUTDIR

# Prepare samples for Bactopia
bactopia prepare \
    --path "$READS_DIR" \
    --species "Mycobacterium bovis" \
    --genome-size 4410000 \
    > "$OUTDIR/Mbovis_samples.txt"

# Run Bactopia pipeline
bactopia \
    --samples "$OUTDIR/Mbovis_samples.txt" \
    --coverage 100 \
    --outdir "$OUTDIR" \
    --max_cpus 8

# Generate summary
bactopia summary \
    --bactopia-path "$OUTDIR"
