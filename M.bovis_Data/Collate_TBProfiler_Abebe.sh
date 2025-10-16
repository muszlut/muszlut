#!/bin/bash
#SBATCH --job-name=Mbovis_CollateTBprofiler
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

# Load Micromamba
module load Micromamba/2.3.0
eval "$(micromamba shell hook --shell bash)"
micromamba activate tbprofiler

# Set directories
OUTDIR="/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples/all_fastqs/TBprofiler_results_conda"
FOFN="/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples/all_fastqs/TBprofiler_results_conda"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Collate TB-Profiler results
#samples just needs to be a list of sample names. No path is required.
tb-profiler collate --samples $FOFN/tbprofiler_results_paths_with_names_absolute.fofn --dir $OUTDIR/*/results --itol
# Deactivate the Micromamba environment
micromamba deactivate



