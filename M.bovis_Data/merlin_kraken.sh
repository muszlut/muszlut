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
OUTDIR="/scratch/ma95362/ETH_bovis_Sequence/Abebe_samples/sub_raw"
DATASETS_DIR="/scratch/ma95362/ETH_bovis_Sequence/Abebe_samples/bactopia-datasets"

# Create output directories if they don't exist
mkdir -p "$OUTDIR"
mkdir -p "$DATASETS_DIR"

# Load Bactopia module
module load Bactopia/3.1.0

# Move to output directory
cd "$OUTDIR"

# (Optional) Download Bactopia datasets if not already done
# You can comment this out if you already have the datasets
bactopia datasets --species "Mycobacterium bovis" --output "$DATASETS_DIR"

# Prepare samples for Bactopia
bactopia prepare \
    --path "$READS_DIR" \
    --species "Mycobacterium bovis" \
    --genome-size 4410000 \
    > "$OUTDIR/Mbovis_samples.txt"

# Run Bactopia with Merlin and Kraken modules enabled
bactopia \
    --samples "$OUTDIR/Mbovis_samples.txt" \
    --datasets "$DATASETS_DIR" \
    --coverage 100 \
    --outdir "$OUTDIR" \
    --max_cpus 8 \
    --wf merlin kraken2

# Generate a summary report
bactopia summary \
    --bactopia-path "$OUTDIR"
