#!/bin/bash
#SBATCH --job-name=CollateTBprofiler          # Job name
#SBATCH --partition=batch                     # Queue/partition
#SBATCH --ntasks=1                            # Number of tasks (1 for single-process)
#SBATCH --cpus-per-task=8                     # Number of CPU cores
#SBATCH --mem=40gb                            # Memory allocation
#SBATCH --time=07-00:00:00                    # Max runtime (7 days)
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out   # STDOUT log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err    # STDERR log
#SBATCH --mail-type=END,FAIL                  # Email notifications
#SBATCH --mail-user=ma95362@uga.edu           # Email recipient

set -euo pipefail  # safer bash settings

# Load Micromamba module
module load Micromamba/2.3.0

# Initialize Micromamba for bash
eval "$(micromamba shell hook --shell bash)"

# Activate TB-Profiler environment
micromamba activate tbprofiler

# Set working directories
OUTDIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/TBprofiler_results_conda"
FOFN="/scratch/ma95362/eth_national_analysis/all_fastq_reads/TBprofiler_results_conda/all_samples.fofn"

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Collate TB-Profiler results
tb-profiler collate --samples all_samples.fofn --dir /scratch/ma95362/eth_national_analysis/all_fastq_reads/TBprofiler_results_conda/*/results --itol


# Deactivate the Micromamba environment
micromamba deactivate