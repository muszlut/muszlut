#!/bin/bash
#SBATCH --job-name=retry_failed_sra
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=06-00:00:00
#SBATCH --output=/scratch/ma95362/EPTB_Hilina/new_project/curl_download_all/log.%j.out
#SBATCH --error=/scratch/ma95362/EPTB_Hilina/new_project/curl_download_all/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

echo "Starting retry at: $(date)"


# Folder with FASTQ files
OUTDIR="/scratch/ma95362/EPTB_Hilina/new_project/curl_download_all"

# List of failed files from previous log:
FAILED_LIST="failed_downloads.txt"

# Make sure the file exists
if [[ ! -f "$FAILED_LIST" ]]; then
    echo "ERROR: $FAILED_LIST file not found!"
    exit 1
fi

cd "$OUTDIR"

while read FILE; do
    if [[ -z "$FILE" ]]; then
        continue
    fi

    # Extract prefix like SRR31229081
    PREFIX=$(echo $FILE | sed 's/_.*//')

    # Construct FTP URL using EBI structure
    SUBDIR=$(echo $PREFIX | sed 's/.\{3\}$//')  # remove last 3 digits
    LAST3=$(echo $PREFIX | tail -c 4)          # last 3 digits

    URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SUBDIR}/${LAST3}/${PREFIX}/${FILE}"

    echo "Retrying: $FILE"
    echo "URL: $URL"

    curl -O "$URL"

    if [[ -f "$FILE" ]]; then
        echo "SUCCESS: $FILE downloaded"
    else
        echo "FAILED AGAIN: $FILE" >> failed_second_round.txt
    fi

done < "$FAILED_LIST"

echo "Completed retry at: $(date)"
