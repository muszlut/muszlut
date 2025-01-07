#!/bin/bash
#SBATCH --job-name=Musse_Spoligotyping                     # Job name
#SBATCH --partition=batch                                  # Partition (queue) name
#SBATCH --ntasks=1                                         # Run on a single CPU
#SBATCH --cpus-per-task=8                                  # Number of cores per task
#SBATCH --mem=16gb                                         # Job memory request
#SBATCH --time=03-00:00:00                                 # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out       # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err        # Standard error log
#SBATCH --mail-type=END,FAIL                               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                        # Where to send mail 

# Variables
# Define the path to your sample sequence reads
SAMPLE_DIR="/work/fdqlab/Ethiopia_wgs_mtb_2024/first_run"  # Replace with the actual path to your sample reads

# Define the path to the output directory for storing results
OUTPUT_DIR="/scratch/ma95362/musse_MGA/fastqs"   # Replace with the actual path for output results

# Define the path to the directory where Bactopia is installed
BACTOPIA_DIR="/apps/mf/eb/all/Bactopia/3.1.0.lua"  # Replace with the actual path to Bactopia installation

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Run Bactopia spoligotyping analysis
bactopia spoligotyping --input $SAMPLE_DIR --output $OUTPUT_DIR --bactopia $BACTOPIA_DIR

# Print a message indicating completion
echo "Spoligotyping analysis complete. Results are in $OUTPUT_DIR"
