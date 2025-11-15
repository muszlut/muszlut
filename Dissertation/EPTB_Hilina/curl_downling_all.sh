#!/bin/bash
#SBATCH --job-name=curl_download_all
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=06-00:00:00
#SBATCH --output=/scratch/ma95362/EPTB_Hilina/new_project/curl_download_all/log.%j.out
#SBATCH --error=/scratch/ma95362/EPTB_Hilina/new_project/curl_download_all/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# -------------------------------
# Step 0: Set Output Directory
# -------------------------------
OUTDIR="/scratch/ma95362/EPTB_Hilina/new_project/curl_download_all"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Log file
LOGFILE="$OUTDIR/download.log"
touch $LOGFILE

echo "Starting Downloads at: $(date)" | tee -a $LOGFILE

# -------------------------------
# Step 1: Define download function
# -------------------------------
download() {
    local URL=$1
    local FILE=$(basename $URL)
    echo "Downloading: $FILE from $URL" | tee -a $LOGFILE
    curl -O -L --retry 5 --retry-delay 10 --continue-at - "$URL"
    if [[ $? -ne 0 ]]; then
        echo "ERROR: Failed to download $FILE" | tee -a $LOGFILE
    else
        echo "SUCCESS: $FILE downloaded" | tee -a $LOGFILE
    fi
}

# -------------------------------
# Step 2: Download all URLs from urls.txt
# -------------------------------
URL_LIST="/scratch/ma95362/EPTB_Hilina/new_project/urls.txt"

if [[ ! -f "$URL_LIST" ]]; then
    echo "ERROR: $URL_LIST not found!" | tee -a $LOGFILE
    exit 1
fi

while read -r URL; do
    # Skip empty lines or lines starting with #
    [[ -z "$URL" || "$URL" == \#* ]] && continue
    download "$URL"
done < "$URL_LIST"

echo "Downloads Completed at: $(date)" | tee -a $LOGFILE
