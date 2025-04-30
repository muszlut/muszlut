#!/bin/bash
#SBATCH --job-name=mtbvartools_pipeline               # Job name
#SBATCH --partition=batch                             # Partition (queue) name
#SBATCH --ntasks=1                                    # Run on a single CPU
#SBATCH --cpus-per-task=8                             # Number of cores per task
#SBATCH --mem=40gb                                    # Job memory request
#SBATCH --time=07-00:00:00                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out  # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err   # Standard error log

#SBATCH --mail-type=END,FAIL                          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                   # Where to send mail



# Load required modules (adjust to your environment)
module purge
module load Python/3.10.4-GCCcore-11.3.0
module load BWA
module load SAMtools
module load GATK
module load TBProfiler
module load FastQC
module load R

# Activate virtual environment if you have one
# source /path/to/your/env/bin/activate

# Optional: activate conda if used
# module load Anaconda3
# conda activate mtbvartools_env

# Define your input/output
# Paths
OUTPUT_DIR="/scratch/ma95362/Sequence/mtbvartools_output"
SCRIPT_PATH="/scratch/ma95362/mtbvartools/scripts/sra_download.py"  # ← Update this path

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

FASTA_REF=/path/to/ref.fasta
GENBANK_REF=/path/to/ref.gbk
FASTQ_PATH=/path/to/sample_1.fastq.gz,/path/to/sample_2.fastq.gz
OUTPUT_NAME=sample001


# Run the script
python3 mtbvartools_pipeline.py \
  --fastq-path "$FASTQ_PATH" \
  --fasta "$FASTA_REF" \
  --genbank "$GENBANK_REF" \
  --output "$OUTPUT_NAME" \
  --dir "$OUTPUT_DIR" \
  --threads 8 \
  --memory 64000m \
  --target-depth 100 \
  --tbprofiler-fastq \
  --overwrite
