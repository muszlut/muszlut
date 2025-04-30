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

# Load correct SPAdes module
module load SPAdes/3.15.5-GCC-11.3.0

# Set working directory
cd /scratch/ma95362/Sequence/Ref_H37Rv/sra_download

# Compress if needed (safe to re-run)
gzip -f sample_ERR2679299_1.fastq
gzip -f sample_ERR2679299_2.fastq

# Run SPAdes
spades.py \
  -1 sample_ERR2679299_1.fastq.gz \
  -2 sample_ERR2679299_2.fastq.gz \
  -o spades_output_ERR2679299 \
  --threads 16 \
  --memory 64

echo "Assembly complete. Output located in spades_output_ERR2679299/contigs.fasta"
