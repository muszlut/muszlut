#!/bin/bash
#SBATCH --job-name=mtbvartools_pipeline               # Job name
#SBATCH --partition=batch                             # Partition (queue) name
#SBATCH --ntasks=1                                    # Run on a single CPU
#SBATCH --cpus-per-task=16                            # Number of cores per task
#SBATCH --mem=64gb                                    # Job memory request
#SBATCH --time=07-00:00:00                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out  # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err   # Standard error log

#SBATCH --mail-type=END,FAIL                          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                   # Where to send mail

# Load required modules
module purge
module load Python/3.10.4-GCCcore-11.3.0
module load BWA
module load SAMtools
module load GATK
module load TBProfiler
module load FastQC
module load R

# Activate virtual environment
source ~/.bashrc
conda activate mtbvartools

# Define paths and inputs
OUTPUT_DIR="/scratch/ma95362/Sequence/mtbvartools_output"
SCRIPT_PATH="/scratch/ma95362/mtbvartools/scripts/sra_variant_pipeline.py"

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

FASTA_REF="/scratch/ma95362/Sequence/Ref_H37Rv/sra_download/spades_output_ERR2679299/contigs.fasta"
GENBANK_REF="/scratch/ma95362/gbk/ncbi_dataset/data/GCF_000195955.2.gbk"
FASTQ_PATH="/scratch/ma95362/Sequence"
OUTPUT_NAME="sample001"

# Run the MTBVarTools pipeline
python3 "$SCRIPT_PATH" \
  --fastq-path "$FASTQ_PATH" \
  --fasta "$FASTA_REF" \
  --genbank "$GENBANK_REF" \
  --output "$OUTPUT_NAME" \
  --dir "$OUTPUT_DIR" \
  --threads 16 \
  --memory 64000m \
  --target-depth 100 \
  --tbprofiler-fastq \
  --overwrite
