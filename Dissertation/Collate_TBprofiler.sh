#!/bin/bash
#SBATCH --job-name=CollateTBprofiler
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40gb
#SBATCH --time=07-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

set -euo pipefail

module load Micromamba/2.3.0
eval "$(micromamba shell hook --shell bash)"
micromamba activate tbprofiler

# Paths
OUTDIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/TBprofiler_results_conda"
FOFN="/scratch/ma95362/eth_national_analysis/all_fastq_reads/all_samples.fofn"

# Ensure FOFN exists
if [ ! -f "$FOFN" ]; then
    echo "FOFN file not found: $FOFN"
    exit 1
fi

# Create output directory
mkdir -p "$OUTDIR"

# Collate TB-Profiler results using the FOFN
tb-profiler collate --samples "$FOFN" --itol

# Deactivate environment
micromamba deactivate
