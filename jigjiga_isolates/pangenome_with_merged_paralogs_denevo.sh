#!/bin/bash
#SBATCH --job-name=pangenome_denevo_mergedparalogs
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120gb
#SBATCH --time=07-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# -------------------------------
# Step 0: Set Output Directory
# -------------------------------

OUTDIR="/scratch/ma95362/publication"
# -------------------------------
# Load modules
# -------------------------------

module load Bactopia/4.0.0-conda
module load Java/17.0.6

# -------------------------------
# Prepare workspace
# -------------------------------

mkdir -p "$OUTDIR"
cd $OUTDIR

# -------------------------------
# Run Bactopia Pangenome
# -------------------------------

bactopia \
    --wf pangenome \
    --bactopia $OUTDIR \
    --panaroo_merge_paralogs