#!/bin/bash
#SBATCH --job-name=TBprofiler_Conda           # Job name
#SBATCH --partition=batch                      # Queue/partition
#SBATCH --ntasks=1                             # Number of tasks (1 for single-process)
#SBATCH --cpus-per-task=8                      # Number of CPU cores
#SBATCH --mem=40gb                             # Memory allocation
#SBATCH --time=07-00:00:00                     # Max runtime (7 days)
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out   # STDOUT log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err    # STDERR log
#SBATCH --mail-type=END,FAIL                   # Email notifications
#SBATCH --mail-user=ma95362@uga.edu           # Email recipient

set -euo pipefail  # safer bash settings

# Load Conda module
module load Mamba/23.1.0-4  # or the available conda/mamba module on your HPC

# Activate TB-Profiler conda environment
# Only need to create it once; afterwards you just activate
CONDA_ENV="$HOME/.conda/envs/tbprofiler"
if [ ! -d "$CONDA_ENV" ]; then
    mamba create -y -n tbprofiler tb-profiler=6.6.5 samtools bcftools -c bioconda
fi
source activate tbprofiler

# Set working directories
FASTQ_DIR="/scratch/ma95362/musse_MGA/merged/first_run_merged/all_reads"
OUTDIR="/scratch/ma95362/my_tbprofiler_results"

mkdir -p "$OUTDIR"

# Change to FASTQ directory
cd "$FASTQ_DIR"

# Loop through all samples
for sample in *._R1.fastq.gz; do
    base=$(basename "$sample" ._R1.fastq.gz)
    SAMPLE_OUT="$OUTDIR/$base"
    mkdir -p "$SAMPLE_OUT"
    
    tb-profiler profile \
        -1 "${base}._R1.fastq.gz" \
        -2 "${base}._R2.fastq.gz" \
        -p "$base" \
        --dir "$SAMPLE_OUT" \
        --threads 8 \
        --txt
done
# Deactivate the conda environment
source deactivate