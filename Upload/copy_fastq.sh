#!/bin/bash
#SBATCH --job-name=copy_clean_reads
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40gb
#SBATCH --time=05-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/copy_clean_reads.%j.out
#SBATCH --error=/scratch/ma95362/scratch/copy_clean_reads.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

module load parallel

SRC=/work/fdqlab/Ethiopia_dataset_763
DEST=/scratch/ma95362/clean_sequences_reads

mkdir -p "$DEST"
cd "$SRC" || { echo "Failed to cd into $SRC"; exit 1; }

# Copy only non-ERR and non-SRR fastq.gz files
find "$SRC" -maxdepth 1 -type f -name "*.fastq.gz" \
  ! -name "ERR*" ! -name "SRR*" \
  | parallel -j 8 rsync -avh --partial {} "$DEST/"
