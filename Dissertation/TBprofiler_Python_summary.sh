#!/bin/bash
#SBATCH --job-name=Summary_TBProfiler        # Job name
#SBATCH --partition=batch                    # Queue/partition
#SBATCH --ntasks=1                           # Single task
#SBATCH --cpus-per-task=8                    # CPU cores
#SBATCH --mem=40gb                           # Memory allocation
#SBATCH --time=02-00:00:00                   # Max runtime (2 days)
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

set -euo pipefail

# Load Micromamba
module load Micromamba/2.3.0

# Set Micromamba root manually (adjust to your actual path if different)
export MAMBA_ROOT_PREFIX=/home/ma95362/micromamba

# Initialize Micromamba for bash
eval "$(micromamba shell hook --shell bash)"

# Activate TB-Profiler environment
micromamba activate tbprofiler

# Define directories
RESULTS_DIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/TBprofiler_results_conda"
OUTPUT_CSV="/scratch/ma95362/eth_national_analysis/all_fastq_reads/TBprofiler_summary.csv"

# Run the Python merge script
python /home/ma95362/muszlut/all_Bactopia_Analysis/python_tbprofiler_custom_script.py \
    --dir "$RESULTS_DIR" \
    --out "$OUTPUT_CSV"

# Deactivate Micromamba
micromamba deactivate
