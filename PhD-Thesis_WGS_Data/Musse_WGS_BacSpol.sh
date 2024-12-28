#!/bin/bash
#SBATCH --job-name=Musse_WGS_Spoligotyping                     # Job name
#SBATCH --partition=batch                                      # Partition (queue) name
#SBATCH --ntasks=1                                             # Run on a single CPU
#SBATCH --cpus-per-task=8                                      # Number of cores per task
#SBATCH --mem=40gb                                             # Job memory request
#SBATCH --time=03-00:00:00                                     # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err            # Standard error log

#SBATCH --mail-type=END,FAIL                                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                            # Email for notifications

# Define directories
READSDIR="/work/fdqlab/Ethiopia_wgs_mtb_2024/first_run"        # Directory containing FASTQ files
DBDIR="/scratch/ma95362/musse_MGA/fastqs/MGA_paired_end_samples/bactopia-runs/bactopia-20241218-162857/merged-results"   # Path to the Bactopia database
OUTDIR="/scratch/ma95362/musse_MGA/fastqs/MGA_paired_end_samples/spoligotyping_results"      # Output directory for Bactopia results

# Ensure output directory exists
mkdir -p "$OUTDIR"

# Load necessary modules
module load Bactopia/3.1.0                                    # Ensure correct module is loaded

# Activate conda environment for Bactopia (if needed)
# module load Miniconda3
# conda activate bactopia-env

# Run Bactopia pipeline with spoligotyping and TB-Profiler enabled
bactopia \
    --path "$READSDIR" \
    --datasets "$DBDIR" \
    --outdir "$OUTDIR" \
    --tools tb-profiler \
    --skip_qc \
    --conda-env "$(conda info --base)/envs/bactopia-env"

# Completion message
echo "Bactopia spoligotyping analysis completed successfully."
