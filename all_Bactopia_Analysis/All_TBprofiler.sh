#!/bin/bash
#SBATCH --job-name=TBprofiler_Run                    # Job name
#SBATCH --partition=batch                            # Queue/partition
#SBATCH --ntasks=1                                  # Number of tasks (1 for single-process)
#SBATCH --cpus-per-task=8                           # Number of CPU cores
#SBATCH --mem=40gb                                  # Memory allocation
#SBATCH --time=07-00:00:00                          # Max runtime (7 days)
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out  # STDOUT log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err   # STDERR log

#SBATCH --mail-type=END,FAIL                         # Email notifications
#SBATCH --mail-user=ma95362@uga.edu                  # Email recipient

set -euo pipefail  # safer bash settings

# Paths
OUTDIR="/scratch/ma95362/my_tbprofiler_results"
FASTQ_DIR="/scratch/ma95362/musse_MGA/merged/first_run_merged/all_reads"

# Create output directory if it does not exist
mkdir -p "$OUTDIR"

# Load TB-Profiler module
module load TB-Profiler/6.6.5

# Change working directory to FASTQ folder
cd "$FASTQ_DIR"

# Loop through all R1 FASTQ files
for sample in *._R1.fastq.gz
do
    base=$(basename "$sample" ._R1.fastq.gz)
    
    tb-profiler profile \
        -1 "${base}._R1.fastq.gz" \
        -2 "${base}._R2.fastq.gz" \
        -p "$base" \
        --dir "${OUTDIR}/${base}" \
        --threads 8 \
        --txt
done
