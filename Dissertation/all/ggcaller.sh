#!/bin/bash
#SBATCH --job-name=ggcaller_all
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/ma95362/ggcaller_out/log.%j.out
#SBATCH --error=/scratch/ma95362/ggcaller_out/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

echo "Activating ggcaller environment..."
module purge
source ~/.bashrc
conda activate ggcaller-env

# Paths
READDIR="/scratch/ma95362/ggcaller_reads"
OUTDIR="/scratch/ma95362/ggcaller_out"
READLIST="${OUTDIR}/reads_list.txt"

mkdir -p "$OUTDIR"

echo "Generating reads list..."
rm -f "$READLIST"

# Auto-generate list in format: R1,R2
for R1 in ${READDIR}/*_R1*.fastq.gz; do
    R2=${R1/_R1/_R2}
    if [[ -f "$R2" ]]; then
        echo "${R1},${R2}" >> "$READLIST"
    else
        echo "WARNING: Missing R2 for $R1"
    fi
done

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
