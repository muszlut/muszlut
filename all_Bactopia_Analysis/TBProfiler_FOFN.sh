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

set -euo pipefail

module load gsl/2.1-8
module load ncurses/6.5
module load TB-Profiler/6.6.5

# Ensure TB-Profiler libraries are used
export LD_LIBRARY_PATH=/apps/eb/TB-Profiler/6.6.5/lib:$LD_LIBRARY_PATH

OUTDIR="/scratch/ma95362/my_tbprofiler_results"
mkdir -p "$OUTDIR"
cd "$OUTDIR"
FOFN="/scratch/ma95362/musse_MGA/all_reads_Bactopia_Analysis/samples_clean.fofn"

while read -r sample read1 read2; do
    echo "Processing $sample..."
    mkdir -p "$OUTDIR/$sample"

    tb-profiler profile \
        -1 "$read1" \
        -2 "$read2" \
        -p "$sample" \
        --dir "$OUTDIR/$sample" \
        --threads 8
done < "$FOFN"