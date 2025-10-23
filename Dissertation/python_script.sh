#!/bin/bash
#SBATCH --job-name=python_script                             # Job name
#SBATCH --partition=batch                                    # Partition (queue) name
#SBATCH --ntasks=1                                           # Run on a single CPU
#SBATCH --cpus-per-task=8                                    # Number of cores per task
#SBATCH --mem=40gb                                           # Job memory request
#SBATCH --time=03:00:00                                      # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out         # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err          # Standard error log

#SBATCH --mail-type=END,FAIL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                          # Where to send mail
# ------------------------------
# 1. Activate Conda environment
# ------------------------------
module load Miniforge3
source activate pyseer-env
# ------------------------------
# 2. Run the filtering and extraction script
# Set working variables
OUTDIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/pangenome_tools_results/bactopia/bactopia-runs/pangenome_of_L4.2.2.2.2/panaroo/filtered_output/pyseer_drugtype_out"
SCRIPT="/home/ma95362/muszlut/Dissertation/filter_extract_DR.py"

# Load necessary modules
#module load Biopython/1.84-foss-2023b
#module load Python/3.11.5-GCCcore-13.2.0
#
# Upgrade pip and install required Python packages in user space
#pip install --user --upgrade pip
#pip install --user tqdm pathogenprofiler
#
# Move to the working directory
cd $OUTDIR

# Run the script
~/.conda/envs/pyseer-env/bin/python $SCRIPT 