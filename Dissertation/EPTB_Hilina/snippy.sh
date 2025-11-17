#!/bin/bash
#SBATCH --job-name=Snippy
#SBATCH --partition=batch 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/ma95362/EPTB_Hilina/new_project_logs/curl_download_all/Bactopia_prepare/log/log.%j.out
#SBATCH --error=/scratch/ma95362/EPTB_Hilina/new_project_logs/curl_download_all/Bactopia_prepare/log/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# -----------------------------
# Load environment modules
# -----------------------------
module load Bactopia/3.2.0-conda
module load Java/17.0.6


#Set output directory variable
OUTDIR="/scratch/ma95362/EPTB_Hilina/new_project_logs/curl_download_all/Bactopia_prepare"
REF="/scratch/ma95362/gbk/ncbi_dataset/data/GCF_000195955.2/genomic.gbk"

# Create output folder if not exist
mkdir -p $OUTDIR
cd $OUTDIR
# Run snippy workflow using all samples
bactopia \
    --wf snippy \
    --reference "$REF" \
    --bactopia "$OUTDIR" \
    --exclude /scratch/ma95362/EPTB_Hilina/new_project_logs/curl_download_all/Bactopia_prepare/bactopia-exclude.tsv