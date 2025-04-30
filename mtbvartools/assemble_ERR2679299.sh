#!/bin/bash
#SBATCH --job-name=assemble_ERR2679299
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=07-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/logs/assemble_ERR2679299_%j.out
#SBATCH --error=/scratch/ma95362/scratch/logs/assemble_ERR2679299_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Stop the script on error
set -e

# Load SPAdes module
module load spades/3.15.5

# Set working directory to where the FASTQ files are located
WORKDIR="/scratch/ma95362/Sequence/Ref_H37Rv/sra_download"
cd "$WORKDIR"

# Optional: compress FASTQ files only if they are not already gzipped
if [ -f sample_ERR2679299_1.fastq ] && [ ! -f sample_ERR2679299_1.fastq.gz ]; then
    gzip -f sample_ERR2679299_1.fastq
fi

if [ -f sample_ERR2679299_2.fastq ] && [ ! -f sample_ERR2679299_2.fastq.gz ]; then
    gzip -f sample_ERR2679299_2.fastq
fi

# Create output directory for SPAdes assembly
OUTDIR="spades_output_ERR2679299"
mkdir -p "$OUTDIR"

# Run SPAdes with paired-end reads
spades.py \
  -1 sample_ERR2679299_1.fastq.gz \
  -2 sample_ERR2679299_2.fastq.gz \
  -o "$OUTDIR" \
  --threads 16 \
  --memory 64

echo "Assembly completed. Resulting contigs are in: $WORKDIR/$OUTDIR/contigs.fasta"
