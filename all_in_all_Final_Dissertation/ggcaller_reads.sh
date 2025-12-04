#!/bin/bash
#SBATCH --job-name=ggcaller_all_samples
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/ma95362/all_in_all_reads/logs/ggcaller_all.out
#SBATCH --error=/scratch/ma95362/all_in_all_reads/logs/ggcaller_all.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# -------------------------------
# Load environment
# -------------------------------
module purge
module load ggCaller/1.4.1

# -------------------------------
# Define paths
# -------------------------------
READS_LIST=/scratch/ma95362/all_in_all_reads/small_reads.txt   # list of all paired-end samples
REFS_LIST=/scratch/ma95362/all_in_all_reads/refs_list.txt      # curated reference genomes
BALROG_DB=/scratch/ma95362/ggcaller_db
OUTDIR=/scratch/ma95362/all_in_all_reads/ggcaller_results

# -------------------------------
# Prepare output directory
# -------------------------------
mkdir -p "$OUTDIR"
cd "$OUTDIR" || { echo "Failed to cd into $OUTDIR"; exit 1; }

# -------------------------------
# Start logging
# -------------------------------
echo "Job started at $(date)"

# -------------------------------
# Run ggCaller
# -------------------------------
ggcaller \
    --reads "${READS_LIST}" \
    --refs "${REFS_LIST}" \
    --annotation ultrasensitive \
    --alignment pan \
    --aligner def \
    --balrog-db "${BALROG_DB}" \
    --threads 32 \
    --save \
    --out "${OUTDIR}"

# -------------------------------
# Check if ggCaller ran successfully
# -------------------------------
if [ $? -ne 0 ]; then
    echo "ggCaller failed!"
    exit 1
fi

# -------------------------------
# End logging
# -------------------------------
echo "Job finished at $(date)"
