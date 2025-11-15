#!/bin/bash
#SBATCH --job-name=curl_download
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=03-00:00:00
#SBATCH --output=/scratch/ma95362/EPTB_Hilina/new_project_logs/log.%j.out
#SBATCH --error=/scratch/ma95362/EPTB_Hilina/new_project_logs/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# -------------------------------
# Step 0: Set Output Directory
# -------------------------------
OUTDIR="/scratch/ma95362/EPTB_Hilina/new_project/curl_download"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "Starting Downloads at: $(date)"

# A helper function for safe downloads
download() {
    local URL=$1
    echo "Downloading: $URL"
    curl -O -L --retry 5 --retry-delay 10 --continue-at - "$URL"
    if [[ $? -ne 0 ]]; then
        echo "ERROR: Failed to download $URL" >&2
    fi
}

# -------------------------------
# SRR31229007
# -------------------------------
download ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR312/007/SRR31229007/SRR31229007.fastq.gz
download ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR312/007/SRR31229007/SRR31229007_1.fastq.gz
download ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR312/007/SRR31229007/SRR31229007_2.fastq.gz

# -------------------------------
# SRR31229016
# -------------------------------
download ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR312/016/SRR31229016/SRR31229016.fastq.gz
download ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR312/016/SRR31229016/SRR31229016_1.fastq.gz
download ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR312/016/SRR31229016/SRR31229016_2.fastq.gz

echo "Downloads Completed at: $(date)"
