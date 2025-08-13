#!/bin/bash
#SBATCH --job-name=TBprofiler_Mamba           # Job name
#SBATCH --partition=batch                     # Queue/partition
#SBATCH --ntasks=1                            # Number of tasks (1 for single-process)
#SBATCH --cpus-per-task=8                     # Number of CPU cores
#SBATCH --mem=40gb                            # Memory allocation
#SBATCH --time=07-00:00:00                    # Max runtime (7 days)
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out   # STDOUT log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err    # STDERR log
#SBATCH --mail-type=END,FAIL                  # Email notifications
#SBATCH --mail-user=ma95362@uga.edu          # Email recipient

set -euo pipefail  # safer bash settings

# Load Micromamba module
module load Micromamba/2.3.0

# Initialize Micromamba for bash
eval "$(micromamba shell hook --shell bash)"

# Activate TB-Profiler environment
micromamba activate tbprofiler

# Set working directories
FASTQ_DIR="/scratch/ma95362/musse_MGA/merged/first_run_merged/all_reads"
OUTDIR="/scratch/ma95362/my_tbprofiler_results"
FOFN="/scratch/ma95362/musse_MGA/all_reads_Bactopia_Analysis/samples_clean.fofn"

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Loop through each line in the FOFN (sample, read1, read2)
while read -r sample read1 read2; do
    echo "Processing $sample..."
    SAMPLE_OUT="$OUTDIR/$sample"
    mkdir -p "$SAMPLE_OUT"

    tb-profiler profile \
        -1 "$read1" \
        -2 "$read2" \
        -p "$sample" \
        --db /home/ma95362/.conda/envs/mtbvartools/share/tbprofiler/tbdb \
        --spoligotype \
        --dir "$SAMPLE_OUT" \
        --threads 8 \
        --txt
done < "$FOFN"

# Deactivate the Micromamba environment
micromamba deactivate
