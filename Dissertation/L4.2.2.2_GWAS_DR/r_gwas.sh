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

# Set output directory
OUTDIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/pangenome_tools_results/bactopia/bactopia-runs/pangenome_of_L4.2.2.2.2/panaroo/filtered_output/pyseer_out"
TIDY="/home/ma95362/muszlut/Dissertation/L4.2.2.2_GWAS_DR/pyseer_analysis.R"

# Create output directory if it doesn't exist
if [ ! -d "$OUTDIR" ]; then
    mkdir -p "$OUTDIR"
fi

cd "$OUTDIR"

# Activate Conda environment
source /apps/conda/etc/profile.d/conda.sh
conda activate r-tidyverse
export R_LIBS_USER=/home/ma95362/.conda/envs/r-tidyverse/lib/R/library
# Load R module
#module load R/4.4.1-gfbf-2023b
# Run the R script
Rscript "$TIDY"
