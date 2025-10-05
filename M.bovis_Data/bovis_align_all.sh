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

# Load conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bovisanalyzer


# Load required modules
module load bwa
module load samtools

# Define directories
FASTQ_DIR="/scratch/ma95362/ETH_bovis_Sequence/bovisanalyzer_output/rasusa"
REF="/scratch/ma95362/ETH_bovis_Sequence/bovisanalyzer_output/bwa/index/bwa/AF2122_97"
BAM_DIR="/scratch/ma95362/ETH_bovis_Sequence/bovisanalyzer_output/samtools"

# Make output directory if not exists
mkdir -p $BAM_DIR

cd $BAM_DIR

# Loop through all forward reads
for R1 in ${FASTQ_DIR}/*_1.fastq.gz; do
    R2=${R1/_1.fastq.gz/_2.fastq.gz}
    SAMPLE=$(basename $R1 _1.fastq.gz)

    # Skip if BAM already exists
    if [[ -f ${BAM_DIR}/${SAMPLE}.sorted.bam ]]; then
        echo "Skipping $SAMPLE (already processed)"
        continue
    fi

    echo "Processing sample: $SAMPLE"

    # Align, convert to BAM, sort, and index
    bwa mem -t 16 $REF $R1 $R2 | \
        samtools view -bS - | \
        samtools sort -o ${BAM_DIR}/${SAMPLE}.sorted.bam

    samtools index ${BAM_DIR}/${SAMPLE}.sorted.bam
done

echo "✅ Alignment complete for all samples."
