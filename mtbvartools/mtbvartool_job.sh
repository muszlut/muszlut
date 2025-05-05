#!/bin/bash
#SBATCH --job-name=mtbvartools_run
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=07-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out  # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err   # Standard error log


# Load required modules
module purge  # Clean environment
module load GCC/12.3.0  # Load compatible GCC version (matches Python/3.11.3-GCCcore-12.3.0 if needed)
# Only load Python if you're NOT using conda. If using conda, skip the Python module line.

# Activate conda environment
source ~/.bashrc  # Ensure conda is available
conda activate mtbvartools

# Move to your working directory (if needed)
cd $SLURM_SUBMIT_DIR  # Automatically go to where you submitted the job

# Run your mtbvartools command
mtbvartools run --input my_samples.csv --outdir output_dir --threads 4
