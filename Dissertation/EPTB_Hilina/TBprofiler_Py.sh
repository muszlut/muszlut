#!/bin/bash
#SBATCH --job-name=Tbprofiler_jason_collate_conversion       # Job name
#SBATCH --partition=batch                                    # Partition (queue) name
#SBATCH --ntasks=1                                           # Run on a single CPU
#SBATCH --cpus-per-task=8                                    # Number of cores per task
#SBATCH --mem=40gb                                           # Job memory request
#SBATCH --time=03:00:00                                      # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out         # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err          # Standard error log

#SBATCH --mail-type=END,FAIL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                          # Where to send mail

# Set working variables
OUTDIR="/scratch/ma95362/EPTB_Hilina/Bactopia_Run/split_reads/TBprofiler_results_conda"
SCRIPT="/home/ma95362/muszlut/Dissertation/EPTB_Hilina/TBprofiler.py"

# Load necessary modules
module load Biopython/1.84-foss-2023b
module load Python/3.11.5-GCCcore-13.2.0
module load TB-Profiler/6.6.5

# Upgrade pip and install required Python packages in user space
pip install --user --upgrade pip
pip install --user tqdm pathogenprofiler

# Move to the working directory
cd $OUTDIR

# Run the script
python $SCRIPT --dir results/
