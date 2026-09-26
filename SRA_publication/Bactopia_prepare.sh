#!/bin/bash
#SBATCH --job-name=Bactopia_prepare_&_run
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=07-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

set -euo pipefail

# Load Mamba
module load Mamba/23.11.0-0

# Activate Bactopia environment
source activate bactopia4

# -------------------------------
# Variables
# -------------------------------
READS_DIR="/scratch/ma95362/clean_sequences_reads"
OUTDIR="/scratch/ma95362/Bio_project_publication"
#REF="/scratch/ma95362/gbk/ncbi_dataset/data/GCF_000195955.2/genomic.gbk"

# -------------------------------
# Create output directory
# -------------------------------
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# -------------------------------
# Prepare samples
# -------------------------------
#bactopia prepare \
#    --path "$READS_DIR" \
#    --species "Mycobacterium tuberculosis" \
#    --genome-size 4410000 \
#    > "$OUTDIR/samples.txt"

## -------------------------------
# Run Bactopia
# -------------------------------
export NXF_OPTS="-Xms2g -Xmx8g"

bactopia \
    --samples "$OUTDIR/samples.txt" \
    --coverage 100 \
    --max_cpus 16 \
    --outdir "$OUTDIR" \
    -resume \
    -process.maxForks 4
` 
#bactopia summary --bactopia-path "$OUTDIR"
#bactopia \
#    --wf snippy \
#    --reference $REF \
#    --bactopia $OUTDIR 
#bactopia \
#    --wf pangenome \
#    --bactopia $OUTDIR