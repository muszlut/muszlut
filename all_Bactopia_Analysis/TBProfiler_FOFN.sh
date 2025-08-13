#!/bin/bash
#SBATCH --job-name=TBprofiler_FOFN           # Job name
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

# Load TB-Profiler module
module load TB-Profiler/6.6.5

# Set directories
FOFN="/scratch/ma95362/musse_MGA/all_reads_Bactopia_Analysis/samples_clean.fofn"
OUTDIR="/scratch/ma95362/my_tbprofiler_results"

mkdir -p "$OUTDIR"

# Change to the output directory
cd "$OUTDIR"

# Loop through samples in FOFN
while read -r read1 read2; do
    base=$(basename "$read1" _R1.fastq.gz)
    SAMPLE_OUT="$OUTDIR/$base"
    mkdir -p "$SAMPLE_OUT"

    tb-profiler profile \
        -1 "$read1" \
        -2 "$read2" \
        -p "$base" \
        --dir "$SAMPLE_OUT" \
        --threads 8 \
        --txt
done < "$FOFN"
