#!/bin/bash
#SBATCH --job-name=ggcaller_all_test
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=10:00:00
#SBATCH --output=/scratch/ma95362/test_ggcaller_out/log.%j.out
#SBATCH --error=/scratch/ma95362/test_ggcaller_out/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

echo "Activating ggcaller environment..."
module purge
source ~/.bashrc
conda activate ggcaller-env

# Paths
#READDIR="/scratch/ma95362/ggcaller_reads"
OUTDIR="/scratch/ma95362/test_ggcaller_out"
READLIST="${OUTDIR}/test_list.txt"

mkdir -p "$OUTDIR"

echo "Generating reads list..."
echo "Reads list generated:"
cat "$READLIST"
echo "Running ggCaller..."
ggcaller \
    --refs "$READLIST" \
    --out "$OUTDIR" \
    --threads 16 \
    --clean-mode strict \
    --annotation fast \
    --identity-cutoff 0.98 \
    --len-diff-cutoff 0.98 \
    --family-threshold 0.7 \
    --core-threshold 0.99

echo "Job completed."
