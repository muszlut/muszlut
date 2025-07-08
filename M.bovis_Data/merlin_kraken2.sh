#!/bin/bash
#SBATCH --job-name=Musse_Mbovis_Bactopia                     # Job name
#SBATCH --partition=batch                                    # Partition (queue) name
#SBATCH --ntasks=1                                           # Run on a single CPU
#SBATCH --cpus-per-task=8                                    # Number of cores per task
#SBATCH --mem=40gb                                           # Job memory request
#SBATCH --time=05-00:00:00                                   # Time limit d-hh:mm:ss
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out         # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err          # Standard error log
#SBATCH --mail-type=END,FAIL                                 # Mail notifications
#SBATCH --mail-user=ma95362@uga.edu                          # Your email address

set -euo pipefail

# === Variables ===
READS_DIR="/scratch/ma95362/ETH_bovis_Sequence/reads/64_SH_Sequence_data/raw/sub_raw"
OUTDIR="/scratch/ma95362/ETH_bovis_Sequence/Abebe_samples/sub_raw"
DATASETS_DIR="/scratch/ma95362/ETH_bovis_Sequence/Abebe_samples/bactopia-datasets"

mkdir -p "$OUTDIR" "$DATASETS_DIR"

# === Load Bactopia ===
module load Bactopia/3.1.0

cd "$OUTDIR"

echo "=== Checking Bactopia datasets ==="
echo "⬇️  (Re)Downloading datasets to $DATASETS_DIR with --force..."
bactopia datasets --species "Mycobacterium bovis" --output "$DATASETS_DIR" --force


echo "=== Preparing samples ==="
bactopia prepare \
    --path "$READS_DIR" \
    --species "Mycobacterium bovis" \
    --genome-size 4410000 \
    > "$OUTDIR/Mbovis_samples.txt"

if [ ! -s "$OUTDIR/Mbovis_samples.txt" ]; then
    echo "�� ERROR: No samples detected in $READS_DIR. Mbovis_samples.txt is empty."
    exit 1
fi

echo "=== Running Bactopia workflow (merlin + kraken2) ==="
bactopia \
    --samples "$OUTDIR/Mbovis_samples.txt" \
    --datasets "$DATASETS_DIR" \
    --coverage 100 \
    --outdir "$OUTDIR" \
    --max_cpus 8 \
    --wf merlin kraken2

echo "=== Generating summary report ==="
bactopia summary \
    --bactopia-path "$OUTDIR"

echo "🎉 === Bactopia pipeline completed successfully ==="

