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

# 🧠 Metadata
echo "Running on node: $(hostname)"
echo "SLURM Job ID: $SLURM_JOB_ID"

# 🧬 Robust Conda activation
eval "$(conda shell.bash hook)"
conda activate r-tidyverse

# 📦 Log R version and packages
Rscript -e 'cat("R version:\n"); print(R.version.string); cat("\nLoaded packages:\n"); print(.packages())'

# 📁 Set working directory and script path
WORKDIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/pangenome_tools_results/bactopia/bactopia-runs/pangenome_of_L4.2.2.2.2/R_analysis_output"
SCRIPT="/home/ma95362/muszlut/Dissertation/L4.2.2.2_GWAS_T3ETH/GWIDE_L4.2.2.2.R"

# 🚀 Run R script
cd "$WORKDIR"
Rscript "$SCRIPT"
