#!/bin/bash
#SBATCH --job-name=bovis_alignment
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=05-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# --- Activate Conda Environment ---
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bovisanalyzer

# --- Load Required Modules ---
module load bwa
module load samtools

# --- Define Directories ---
FASTQ_DIR="/scratch/ma95362/ETH_bovis_Sequence/bovisanalyzer_output/rasusa"
REF="/scratch/ma95362/ETH_bovis_Sequence/bovisanalyzer_output/bwa/index/bwa/AF2122_97"
BAM_DIR="/scratch/ma95362/ETH_bovis_Sequence/bovisanalyzer_output/samtools"

# --- Prepare Output Directory ---
mkdir -p $BAM_DIR
cd $BAM_DIR

# --- Alignment Loop ---
for R1 in ${FASTQ_DIR}/*_1.fastq.gz; do
    R2=${R1/_1.fastq.gz/_2.fastq.gz}
    SAMPLE=$(basename $R1 _1.fastq.gz)

    if [[ ! -f "$R2" ]]; then
        echo "❌ Missing pair for $R1 — skipping"
        continue
    fi

    # Skip already processed samples
    if [[ -f "${BAM_DIR}/${SAMPLE}.sorted.bam" ]]; then
        echo "⚠️  Skipping $SAMPLE (already processed)"
        continue
    fi

    echo "🔹 Processing sample: $SAMPLE"

    # Align, convert to BAM, sort, and index
    bwa mem -t 8 $REF $R1 $R2 | \
        samtools view -bS - | \
        samtools sort -@ 8 -o ${BAM_DIR}/${SAMPLE}.sorted.bam -

    samtools index ${BAM_DIR}/${SAMPLE}.sorted.bam
done

echo "✅ Alignment complete for all samples."
