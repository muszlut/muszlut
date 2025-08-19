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
#SBATCH --mail-user=ma95362@uga.edu           # Email recipient

set -euo pipefail  # safer bash settings

# Load Micromamba module
module load Micromamba/2.3.0

# Initialize Micromamba for bash
eval "$(micromamba shell hook --shell bash)"

# Activate TB-Profiler environment
micromamba activate tbprofiler

# Set working directories
OUTDIR="/scratch/ma95362/my_tbprofiler_results"
FOFN="/scratch/ma95362/my_tbprofiler_results/results.fofn"

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Collate TB-Profiler results
tb-profiler collate --samples samples_names.txt --dir /scratch/ma95362/my_tbprofiler_results/*/results --itol


# Deactivate the Micromamba environment
micromamba deactivate
