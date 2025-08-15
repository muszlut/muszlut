#!/bin/bash
#SBATCH --job-name=merge_results                        # Job name
#SBATCH --partition=batch                               # Partition (queue) name
#SBATCH --ntasks=1                                      # Run on a single CPU
#SBATCH --cpus-per-task=8                               # Number of cores per task
#SBATCH --mem=40gb                                      # Job memory request
#SBATCH --time=03:00:00                                 # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out    # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err     # Standard error log

#SBATCH --mail-type=END,FAIL                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                     # Where to send mail	

# Directories
OUTDIR="/scratch/ma95362/my_tbprofiler_results"
SCRIPT="/home/ma95362/muszlut/all_Bactopia_Analysis"

# Load TB-Profiler module (includes pathogenprofiler)
module purge
module load TB-Profiler/6.6.5
#Move to results directory
# Move to results directory
cd $OUTDIR

# Run your Python script
python $SCRIPT/python_tbprofiler_custom_script.py \
    --out merged_results.csv \
    --dir $OUTDIR
