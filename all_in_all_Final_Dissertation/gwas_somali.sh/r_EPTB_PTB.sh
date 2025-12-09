#!/bin/bash
#SBATCH --job-name=R_eptb_ptb
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/ma95362/scratch/R_eptb_ptb_%j.out
#SBATCH --error=/scratch/ma95362/scratch/R_eptb_ptb_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# ------------------------------
# Load Conda and enable "conda activate"
# ------------------------------
module load Miniforge3/24.11.3-0
source $CONDA_PREFIX/etc/profile.d/conda.sh

# ------------------------------
# Activate R environment
# ------------------------------
conda activate r-tidyverse452

# ------------------------------
# Ensure the correct R library path
# ------------------------------
export R_LIBS_USER=$CONDA_PREFIX/lib/R/library

# ------------------------------
# Set directories
# ------------------------------
OUTDIR="/scratch/ma95362/all_in_all_reads/bactopia_prepare/bactopia-runs/pangenome-20251128-070449/panaroo/pyseer_output"
TIDY="/home/ma95362/muszlut/all_in_all_Final_Dissertation/gwas_EPTBvsPTB/pyseer_anlaysis_EPTB_PTB.R"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# ------------------------------
# Run R script
# ------------------------------
Rscript "$TIDY"
