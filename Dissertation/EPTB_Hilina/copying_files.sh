#!/bin/bash
#SBATCH --job-name=copy_SRR3122_dirs
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120gb
#SBATCH --time=03-00:00:00
#SBATCH --output=/scratch/ma95362/EPTB_Hilina/new_project_logs/curl_download_all/Bactopia_prepare/log.%j.out
#SBATCH --error=/scratch/ma95362/EPTB_Hilina/new_project_logs/curl_download_all/Bactopia_prepare/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu


# Define source and destination directories
SRC_DIR="/scratch/ma95362/EPTB_Hilina/new_project_logs/curl_download_all/Bactopia_prepare"
#DEST_DIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples"
DEST_DIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/pangenome_tools_results/bactopia"

# List of directories to exclude
EXCLUDE=("SRR31229027" "SRR31228991" "SRR31229028" "SRR31229109" "SRR31229077" "SRR31229064")
cd $SRC_DIR
echo "Copying SRR3122* directories except excluded ones..."

for dir in "$SRC_DIR"/SRR3122*/; do
    [ -d "$dir" ] || continue
    basename=$(basename "$dir")
    
    # Check if basename is in the exclude list
    skip=false
    for ex in "${EXCLUDE[@]}"; do
        if [[ "$basename" == "$ex" ]]; then
            skip=true
            break
        fi
    done

    if [ "$skip" = false ]; then
        echo "Copying $basename..."
        cp -a "$dir" "$DEST_DIR"/
    else
        echo "Skipping $basename"
    fi
done

echo "Selective copy completed."


