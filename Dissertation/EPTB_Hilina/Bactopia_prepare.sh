#!/bin/bash
#SBATCH --job-name=Bactopia_Mtb
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120gb
#SBATCH --time=03-00:00:00
#SBATCH --output=/scratch/ma95362/EPTB_Hilina/Newe/log.%j.out
#SBATCH --error=/scratch/ma95362/EPTB_Hilina/Newe/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu
# -------------------------------
# ✅ Step 0: (Optional) Download accessions list
# -------------------------------
#To download the accessions list first, uncomment the lines below
#module load Bactopia/3.2.0-conda
#module load EDirect/20.5.20231006-GCCcore-12.3.0
#bactopia search \
#    --query PRJNA1174701 \
#    --force \
#    --genome-size 4400000 \
#    --min-coverage 20
# -------------------------------
# ✅ Step 1: Activate Conda first
# -------------------------------
# Load miniconda module if needed on Sapelo2
module load Miniconda3/22.11.1-1

# Initialize Conda for non-interactive shell
source ~/.bashrc

# Activate your Bactopia Conda environment
conda activate bactopia

# (Alternative — if mamba was installed)
# mamba activate bactopia

# -------------------------------
# ✅ Step 2: Define directories
# -------------------------------
#OUTDIR="/scratch/ma95362/EPTB_Hilina/Newe" This was previous OUTDIR definition until we do the summary step (Step4)
OUTDIR="/scratch/ma95362/EPTB_Hilina/Newe/Bactopia_Run"
# Create output directory
mkdir -p $OUTDIR
cd $OUTDIR
# -------------------------------
# ✅ Step 3: Run Bactopia
# -------------------------------
# Option 1 — run accessions list
#bactopia \
#    --accessions bactopia-accessions.txt \
#    --outdir $OUTDIR/Bactopia_Run \
#    --species "Mycobacterium tuberculosis" \
#    --max_cpus 16
# -------------------------------
# ✅ Step 4: Generate summary
# -------------------------------
bactopia summary \
    --bactopia-path "$OUTDIR"