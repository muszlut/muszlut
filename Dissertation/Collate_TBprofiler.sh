#!/bin/bash
#SBATCH --job-name=CollateTBprofiler          # Job name
#SBATCH --partition=batch                     # Queue/partition
#SBATCH --ntasks=1                            # Single process
#SBATCH --cpus-per-task=8                     # Number of CPU cores
#SBATCH --mem=40gb                            # Memory allocation
#SBATCH --time=07-00:00:00                    # Max runtime (7 days)
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out   # STDOUT log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err    # STDERR log
#SBATCH --mail-type=END,FAIL                  # Email notifications
#SBATCH --mail-user=ma95362@uga.edu           # Email recipient

set -euo pipefail  # safer bash

# Load Micromamba
module load Micromamba/2.3.0
eval "$(micromamba shell hook --shell bash)"
micromamba activate tbprofiler

# Set directories
RESULTS_DIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/TBprofiler_results_conda"
FOFN="/scratch/ma95362/eth_national_analysis/all_fastq_reads/all_samples.fofn"
OUTDIR="${RESULTS_DIR}/collated_results"

mkdir -p "$OUTDIR"

# Generate FOFN with full paths to JSON files
find "$RESULTS_DIR" -type f -name "*.results.json" | sort > "$FOFN"

# Collate TB-Profiler results
tb-profiler collate --samples "$FOFN" --outdir "$OUTDIR" --itol

# Deactivate environment
micromamba deactivate
