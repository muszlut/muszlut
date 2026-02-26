#!/bin/bash
#SBATCH --job-name=TBprofiler_Condas          # Job name
#SBATCH --partition=batch                     # Queue/partition
#SBATCH --ntasks=1                            # Number of tasks (single-process)
#SBATCH --cpus-per-task=16                    # Number of CPU cores
#SBATCH --mem=64gb                            # Memory allocation
#SBATCH --time=07-00:00:00                    # Max runtime (7 days)
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out   # STDOUT log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err    # STDERR log
#SBATCH --mail-type=END,FAIL                  # Email notifications
#SBATCH --mail-user=ma95362@uga.edu          # Email recipient

set -euo pipefail  # safer bash settings

# Load Conda/Miniforge module
module load Miniforge3

# Initialize Conda for bash
eval "$(conda shell.bash hook)"

# Activate TB-Profiler environment
conda activate tbprofiler_env

# Set working directories
READS_DIR="/scratch/ma95362/clean_sequences_reads"
OUTDIR="/scratch/ma95362/publication/tbprofiler_results"
DB="/home/ma95362/.conda/envs/tbprofiler_env/share/tbprofiler/tbdb"

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

cd "$READS_DIR"

# Loop through all R1 files and run tb-profiler
for R1 in *_R1.fastq.gz
do
    SAMPLE=$(basename "$R1" _R1.fastq.gz)
    R2="${SAMPLE}_R2.fastq.gz"

    echo "Processing $SAMPLE"

    tb-profiler profile \
        --read1 "$R1" \
        --read2 "$R2" \
        --threads 16 \
        --spoligotype \
        --db "$DB" \
        --prefix "$SAMPLE" \
        --txt \
        --dir "$OUTDIR"
done