#!/bin/bash
#SBATCH --job-name=TBprofiler_module_load
#SBATCH --partition=highmem_p           # high memory partition for clustering
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32             # parallelize clustering
#SBATCH --mem=240G                     # enough memory for 1398 genomes
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out   # STDOUT log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err    # STDERR log
#SBATCH --mail-type=END,FAIL                  # Email notifications
#SBATCH --mail-user=ma95362@uga.edu          # Email recipient

set -euo pipefail

module load TB-Profiler/6.6.5

# Ensure TB-Profiler libraries are used
export LD_LIBRARY_PATH=/apps/eb/TB-Profiler/6.6.5/lib:$LD_LIBRARY_PATH

OUTDIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/TBprofiler_results_module_load"
mkdir -p "$OUTDIR"
cd "$OUTDIR"
FOFN="/scratch/ma95362/eth_national_analysis/all_fastq_reads/TBprofiler_FOFN_samples.fofn"

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