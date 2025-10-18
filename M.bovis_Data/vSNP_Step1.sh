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

# ─────────────────────────────────────────────
# Load and activate environment
# ─────────────────────────────────────────────
source ~/.bashrc
conda activate vsnp3

# ─────────────────────────────────────────────
# Define input and output directories
# ─────────────────────────────────────────────
FASTQ_DIR="/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples/all_fastqs"
OUT_BASE="$FASTQ_DIR/vSNP_output/step1_output"
EXCLUDE_FILE="$FASTQ_DIR/exclude_list.txt"

mkdir -p "$OUT_BASE"

# ─────────────────────────────────────────────
# Loop through each R1 FASTQ pair
# ─────────────────────────────────────────────
for R1 in ${FASTQ_DIR}/*_R1.fastq.gz; do
    SAMPLE=$(basename "$R1" _R1.fastq.gz)
    R2=${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz
    SAMPLE_OUT="$OUT_BASE/${SAMPLE}"

    # Skip excluded samples
    if grep -q "^${SAMPLE}$" "$EXCLUDE_FILE"; then
        echo "⏩ Skipping excluded sample: $SAMPLE"
        continue
    fi

    echo "▶️ Processing sample: $SAMPLE"
    mkdir -p "$SAMPLE_OUT"
    cd "$SAMPLE_OUT"

    # Clean any leftover partial outputs
    rm -rf sourmash alignment_* unmapped_reads *.bam *.vcf *.csv

    # Run vSNP3 Step 1
    vsnp3_step1.py \
        --FASTQ_R1 "$R1" \
        --FASTQ_R2 "$R2" \
        --output_dir "$SAMPLE_OUT" \
        --spoligo

    echo "✅ Completed sample: $SAMPLE"
    echo "───────────────────────────────────────"
done

echo "🎯 All vSNP Step 1 jobs finished successfully."

