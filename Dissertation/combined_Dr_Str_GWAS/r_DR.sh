#!/bin/bash
#SBATCH --job-name=R_script
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/ma95362/scratch/R_eptb_ptb_%j.out
#SBATCH --error=/scratch/ma95362/scratch/R_eptb_ptb_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Set output directory
OUTDIR="c"
TIDY="/home/ma95362/muszlut/Dissertation/combined_Dr_Str_GWAS/Rscript_for_Lineage_SIT.R"

# Create output directory if it doesn't exist
if [ ! -d "$OUTDIR" ]; then
    mkdir -p "$OUTDIR"
fi

cd "$OUTDIR"

# Activate Conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate r-tidyverse
export R_LIBS_USER=/home/ma95362/.conda/envs/r-tidyve
rse/lib/R/library

# Run the R script
Rscript "$TIDY" > R_output.log 2>&1
