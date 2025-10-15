#!/bin/bash
#SBATCH --job-name=vSNP_step1
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Load and activate environment
source ~/.bashrc
conda activate vsnp3

# Define directories
FASTQ_DIR="/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples/all_fastqs"
OUT_DIR="$FASTQ_DIR/vSNP_output/step1_output"

mkdir -p "$OUT_DIR"

# Load excluded sample names into memory
EXCLUDE_FILE="/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples/all_fastqs/exclude_list.txt"

# Loop through all paired FASTQs
for R1 in ${FASTQ_DIR}/*_R1.fastq.gz; do
    SAMPLE=$(basename "$R1" _R1.fastq.gz)

    # Skip excluded samples
    if grep -q "^${SAMPLE}$" "$EXCLUDE_FILE"; then
        echo "⏩ Skipping excluded sample: $SAMPLE"
        continue
    fi

    R2=${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz
    echo "Processing sample: $SAMPLE"

    vsnp3_step1.py \
        --FASTQ_R1 "$R1" \
        --FASTQ_R2 "$R2" \
        --output_dir "$OUT_DIR" \
        --spoligo
done


echo "✅ All vSNP Step1 jobs finished successfully."
