#!/bin/bash
#SBATCH --job-name=Bactopia_Newproject_Prep
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120gb
#SBATCH --time=03-00:00:00
#SBATCH --output=/scratch/ma95362/EPTB_Hilina/new_project_logs/curl_download_all/log.%j.out
#SBATCH --error=/scratch/ma95362/EPTB_Hilina/new_project_logs/curl_download_all/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# -------------------------------
# Step 0: Set Output Directory
# -------------------------------
OUTDIR="/scratch/ma95362/EPTB_Hilina/new_project_logs/curl_download_all/Bactopia_prepare"
mkdir -p $OUTDIR
cd $OUTDIR
module load Bactopia/3.2.0-conda
bactopia prepare \
    --path /scratch/ma95362/EPTB_Hilina/new_project/curl_download_all \
    --species "Mycobacterium tuberculosis" \
    --genome-size 4400000
    > samples.txt