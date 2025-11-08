#!/bin/bash
#SBATCH --job-name=convert_gzip_pairs
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120GB
#SBATCH --time=05-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/stream_pairs.%j.out
#SBATCH --error=/scratch/ma95362/scratch/stream_pairs.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

module load sratoolkit/3.0.0
module load parallel
module load pigz

DEST=/scratch/ma95362/EPTB_Hilina
cd $DEST

convert_sra() {
    sra=$1
    base=$(basename $sra .sra)

    # Skip if gzipped files already exist
    if [[ -f "${DEST}/${base}_R1.fastq.gz" && -f "${DEST}/${base}_R2.fastq.gz" ]]; then
        echo "Skipping ${base}, already processed."
        return
    fi

    echo "Processing $base..."

    # Convert SRA → FASTQ (temp files go into working dir)
    fasterq-dump --split-files --threads 4 -O . "$sra"

    # Rename to R1/R2
    mv "${base}_1.fastq" "${base}_R1.fastq"
    mv "${base}_2.fastq" "${base}_R2.fastq"

    # Compress with pigz
    pigz -p 4 "${base}_R1.fastq"
    pigz -p 4 "${base}_R2.fastq"

    # Move compressed files into DEST
    mv "${base}_R1.fastq.gz" "$DEST/"
    mv "${base}_R2.fastq.gz" "$DEST/"
}

export -f convert_sra
export DEST

# Run in parallel on remaining .sra files
find $DEST -name "*.sra" | parallel -j 4 convert_sra {}

echo "✅ All SRR files converted into ${DEST}/SRRXXXXXX_R1.fastq.gz and ${DEST}/SRRXXXXXX_R2.fastq.gz"