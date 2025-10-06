#!/bin/bash
#SBATCH --job-name=TBprofiler_conda
#SBATCH --partition=highmem_p           # high memory partition for clustering
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32             # parallelize clustering
#SBATCH --mem=240G                     # enough memory for 1398 genomes
#SBATCH --time=7-00:00:00
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
FASTQ_DIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads"
OUTDIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/TBprofiler_results_conda"
FOFN="/scratch/ma95362/eth_national_analysis/all_fastq_reads/TBprofiler_FOFN_samples.fofn"

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"
cd "$OUTDIR"
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