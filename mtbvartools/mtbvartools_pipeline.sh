#!/bin/bash 
#SBATCH --job-name=mtbvartools_pipeline               # Job name
#SBATCH --partition=batch                             # Partition (queue) name
#SBATCH --ntasks=1                                    # Run on a single CPU
#SBATCH --cpus-per-task=16                            # Number of cores per task
#SBATCH --mem=64gb                                    # Job memory request
#SBATCH --time=07-00:00:00                            # Time limit (7 days)
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out  # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err   # Standard error log

#SBATCH --mail-type=END,FAIL                          # Mail events
#SBATCH --mail-user=ma95362@uga.edu                   # Email address for notifications

# Clean environment and load required modules
module purge
export LMOD_IGNORE_CACHE=yes
module load GCCcore/13.3.0
module load Python/3.11.3-GCCcore-13.3.0
module load BWA
module load SAMtools
module load GATK
module load TBProfiler
module load FastQC
module load R

# Activate conda environment
source ~/.bashrc
source activate mtbvartools

# Add custom script path
export PATH="/scratch/ma95362/mtbvartools/scripts:$PATH"

# Define variables
OUTPUT_DIR="/scratch/ma95362/Sequence/mtbvartools_output"
SCRIPT_PATH="/scratch/ma95362/mtbvartools/scripts/sra_variant_pipeline.py"
FASTA_REF="/scratch/ma95362/Sequence/Ref_H37Rv/sra_download/spades_output_ERR2679299/contigs.fasta"
FASTQ_PATH="/scratch/ma95362/Sequence"
SAMPLE_LIST="/scratch/ma95362/Sequence/samples.txt"
LOG_DIR="/scratch/ma95362/Sequence/logs"

# Create directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"
cd "$OUTPUT_DIR"

# Check if sample list exists
if [[ ! -f "$SAMPLE_LIST" ]]; then
  echo "❌ Sample list not found: $SAMPLE_LIST"
  exit 1
fi

# Run pipeline for each sample
while read SAMPLE; do
  if [[ -f ${FASTQ_PATH}/${SAMPLE}_R1.fastq && -f ${FASTQ_PATH}/${SAMPLE}_R2.fastq ]]; then
    FQ1="${SAMPLE}_R1.fastq"
    FQ2="${SAMPLE}_R2.fastq"
  elif [[ -f ${FASTQ_PATH}/${SAMPLE}_1.fastq && -f ${FASTQ_PATH}/${SAMPLE}_2.fastq ]]; then
    FQ1="${SAMPLE}_1.fastq"
    FQ2="${SAMPLE}_2.fastq"
  else
    echo "❌ Skipping ${SAMPLE}: FASTQ pair not found"
    continue
  fi

  echo "✅ Running pipeline for $SAMPLE..."
  python "$SCRIPT_PATH" \
    --fastq-path "${FASTQ_PATH}/${FQ1},${FASTQ_PATH}/${FQ2}" \
    --reference "$FASTA_REF" \
    --output-dir "${OUTPUT_DIR}/${SAMPLE}" \
    --threads 8 \
    --memory 64000m \
    --target-depth 100 \
    --tbprofiler-fastq \
    --overwrite \
    > "${LOG_DIR}/${SAMPLE}.out" 2> "${LOG_DIR}/${SAMPLE}.err"

done < "$SAMPLE_LIST"
