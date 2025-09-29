#!/bin/bash
#SBATCH --job-name=copy_bactopia
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/ma95362/eth_national_analysis/logs/copy_bactopia_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Load any necessary modules (optional)
# module load rsync

# Source and destination directories
SRC="/scratch/ma95362/eth_national_analysis/bactopia_results"
DEST="/scratch/ma95362/eth_national_analysis/Bactopia_runs"

# Create destination directory if it doesn't exist
mkdir -p "$DEST"

cd "$SRC" || { echo "Failed to cd into $SRC"; exit 1; } 
# Loop over each top-level sample directory
for sample_dir in "$SRC"/*; do
    if [ -d "$sample_dir" ]; then
        sample_name=$(basename "$sample_dir")
        target_subdir="$sample_dir/$sample_name"

        if [ -d "$target_subdir" ]; then
            mkdir -p "$DEST/$sample_name"
            cp -r "$target_subdir"/* "$DEST/$sample_name/"
            echo "Copied $target_subdir to $DEST/$sample_name/"
        else
            echo "No matching subdirectory for $sample_name, skipping."
        fi
    fi
done

echo "Done copying all matching subdirectories."
