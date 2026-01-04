#!/bin/bash
#SBATCH --job-name=prep_shuffled_fasta
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/ma95362/scratch/ncbi_upload.%j.out
#SBATCH --error=/scratch/ma95362/scratch/ncbi_upload.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# -------------------------
# Paths
# -------------------------
FASTQ_DIR=/scratch/ma95362/clean_sequences_reads
OUT_DIR=/scratch/ma95362/CRISPRbuilder-TB/sequences

# -------------------------
# Go to FASTQ directory
# -------------------------
cd $FASTQ_DIR || exit 1

echo "Starting FASTQ → shuffled FASTA preparation"
echo "Input directory: $FASTQ_DIR"
echo "Output directory: $OUT_DIR"

# -------------------------
# Main loop
# -------------------------
for r1 in *_R1.fastq.gz; do
    sample=$(basename "$r1" _R1.fastq.gz)

    echo "--------------------------------------"
    echo "Processing sample: $sample"
    echo "--------------------------------------"

    # Create isolate directory
    mkdir -p $OUT_DIR/$sample

    # Copy paired-end reads
    cp ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz \
       $OUT_DIR/$sample/

    cd $OUT_DIR/$sample || exit 1

    # Unzip FASTQ files
    gunzip -f ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz

    # Convert FASTQ → shuffled FASTA
    awk '{if(NR%4==1){printf(">%s\n",substr($0,2))} else if(NR%4==2){print}}' \
        ${sample}_R1.fastq ${sample}_R2.fastq \
        > ${sample}_shuffled.fasta

    cd $FASTQ_DIR || exit 1
done

echo "All samples processed successfully."
