#!/bin/bash
#SBATCH --job-name=bp_700_mtbc
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120gb
#SBATCH --time=07-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

set -euo pipefail

module load Bactopia/3.2.0-conda
module load Java/17.0.6

OUTDIR="/scratch/ma95362/eth_3rd_national_dataset"
SEARCHDIR="$OUTDIR/search"
RUNNAME="ena_prjna1104194_700"

mkdir -p "$SEARCHDIR"
cd "$OUTDIR"

# Step 1: Search
#bactopia search \
#    --query PRJNA1104194 \
#    --outdir "$SEARCHDIR"

# Step 2: Run
#bactopia \
#    --accessions "$SEARCHDIR/bactopia-accessions.txt" \
#    --coverage 30 \
#    --outdir "$OUTDIR/$RUNNAME" \
#    --max_cpus $SLURM_CPUS_PER_TASK

# Step 3: summary
bactopia summary \
    --bactopia-path "$OUTDIR"