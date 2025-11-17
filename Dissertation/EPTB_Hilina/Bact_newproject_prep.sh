#!/bin/bash
#SBATCH --job-name=Bactopia_Newproject_Prep
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120gb
#SBATCH --time=03-00:00:00
#SBATCH --output=/scratch/ma95362/EPTB_Hilina/new_project_logs/curl_download_all/Bactopia_prepare/log.%j.out
#SBATCH --error=q/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# -------------------------------
# Step 0: Set Output Directory
# -------------------------------
# Variables
#READS_DIR="/scratch/ma95362/EPTB_Hilina/new_project/curl_download_all"
OUTDIR="/scratch/ma95362/EPTB_Hilina/new_project_logs/curl_download_all/Bactopia_prepare"
# Load Bactopia module
module load Bactopia/3.2.0-conda
# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"
# Move to working directory
cd $OUTDIR
# Prepare samples for Bactopia
#bactopia prepare \
#    --path "$READS_DIR" \
#    --species "Mycobacterium tuberculosis" \
#    --genome-size 4410000 \
#    > "$OUTDIR/Hilina_samples.txt"
# Run Bactopia pipeline
#bactopia \
#    --samples "$OUTDIR/Hilina_samples.txt" \
#    --coverage 100 \
#    --outdir "$OUTDIR" \
#    --max_cpus 8
# Generate summary
#bactopia summary \
#    --bactopia-path "$OUTDIR"
bactopia \
    --wf tbprofiler \
    --bactopia $OUTDIR \
    --max_cpus 32