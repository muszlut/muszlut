#!/bin/bash
#SBATCH --job-name=mustutorial         # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=1gb                     # Job memory request
#SBATCH --time=02:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out   # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err    # Standard error log

#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu   # Where to send mail	

# Set output directory variable
OUTDIR="/scratch/ma95362/reference"

# Create the output directory if it doesn't exist
if [ ! -d "$OUTDIR" ]; then 
    mkdir -p "$OUTDIR"
fi

# Load required module
module load NCBI-Datasets-CLI/16.4.4

# Move to working directory
cd "$OUTDIR"

# Download genome dataset with FASTA (genomic) only
datasets download genome accession GCF_000195955.2 --include genome

# Unzip the downloaded dataset
unzip ncbi_dataset.zip