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

echo "Loading Conda module and activating ggcaller environment..."
module purge
module load anaconda3   # Load cluster-provided Conda
conda activate ggcaller-env

# Set TMPDIR to avoid /tmp space issues
export TMPDIR=/scratch/$USER/tmp
mkdir -p "$TMPDIR"

# Paths
OUTDIR="/scratch/ma95362/test_ggcaller_out"
READLIST="${OUTDIR}/test_list.txt"

mkdir -p "$OUTDIR"

# Optional: regenerate reads list automatically
# Find all .fna.gz in /main/assembler directories of your samples
# and save to READLIST
find /scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples \
    -type f -path "*/main/assembler/*.fna.gz" > "$READLIST"

echo "Reads list generated:"
cat "$READLIST"

echo "Running ggCaller..."
ggcaller \
    --reads "$READLIST" \
    --out "$OUTDIR" \
    --threads 16 \
    --clean-mode strict \
    --annotation fast \
    --identity-cutoff 0.98 \
    --len-diff-cutoff 0.98 \
    --family-threshold 0.7 \
    --core-threshold 0.99

echo "Job completed."
