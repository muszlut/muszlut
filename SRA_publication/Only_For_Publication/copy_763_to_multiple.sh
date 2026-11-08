```bash
#!/bin/bash
#SBATCH --job-name=copy_SRR_ERR
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=07-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

set -euo pipefail

# Source and destination
SRC="/scratch/ma95362/Eth_dataset_763/bactopia_run/results"
DEST="/scratch/ma95362/Bio_project_publication/Mtb_Multiple_BioProjects/results_singularity"

# Create destination
mkdir -p "$DEST"

echo "=========================================="
echo "Starting SRR and ERR directory transfer"
echo "Date: $(date)"
echo "Source: $SRC"
echo "Destination: $DEST"
echo "=========================================="

# Identify SRR and ERR directories
shopt -s nullglob

DIRS=("$SRC"/SRR*/ "$SRC"/ERR*/)

SOURCE_COUNT=${#DIRS[@]}

echo "Total sample directories found: $SOURCE_COUNT"

# Copy directories using rsync
for DIR in "${DIRS[@]}"; do

    echo "------------------------------------------"
    echo "Copying: $(basename "$DIR")"

    rsync -a --info=progress2 "$DIR" "$DEST/"

done

echo "=========================================="
echo "Transfer completed: $(date)"
echo "=========================================="

# Count source directories
SOURCE_SRR=$(find "$SRC" -mindepth 1 -maxdepth 1 -type d -name "SRR*" | wc -l)
SOURCE_ERR=$(find "$SRC" -mindepth 1 -maxdepth 1 -type d -name "ERR*" | wc -l)

# Count destination directories
DEST_SRR=$(find "$DEST" -mindepth 1 -maxdepth 1 -type d -name "SRR*" | wc -l)
DEST_ERR=$(find "$DEST" -mindepth 1 -maxdepth 1 -type d -name "ERR*" | wc -l)

echo "========== TRANSFER SUMMARY =============="
echo "Source SRR directories: $SOURCE_SRR"
echo "Destination SRR directories: $DEST_SRR"
echo "Source ERR directories: $SOURCE_ERR"
echo "Destination ERR directories: $DEST_ERR"
echo "=========================================="

if [ "$SOURCE_SRR" -eq "$DEST_SRR" ] && \
   [ "$SOURCE_ERR" -eq "$DEST_ERR" ]; then
    echo "SUCCESS: Directory counts match."
else
    echo "WARNING: Directory counts do not match."
fi

echo "Job finished: $(date)"
```