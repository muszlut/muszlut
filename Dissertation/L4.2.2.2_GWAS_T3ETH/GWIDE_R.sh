#!/bin/bash
#SBATCH --job-name=R_analysis_GWIDE
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/ma95362/scratch/R_analysis_%j.out
#SBATCH --error=/scratch/ma95362/scratch/R_analysis_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Initialize Conda for non-interactive shell
source ~/.bashrc
# Activate Conda
conda activate r-tidyverse
export R_LIBS_USER=/home/ma95362/.conda/envs/r-tidyverse/lib/R/library
# Set working directory
WORKDIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/pangenome_tools_results/bactopia/bactopia-runs/pangenome_of_L4.2.2.2.2/R_analysis_output"
SCRIPT="/home/ma95362/muszlut/Dissertation/L4.2.2.2_GWAS_T3ETH/GWIDE_L4.2.2.2.R"
# Run R script
cd "$WORKDIR"
Rscript "$SCRIPT"
