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

# Ensure the output directory exists
if [ ! -d "$OUTDIR" ]; then
    mkdir -p "$OUTDIR"
fi

# Load Singularity module (if required by the HPC system)
module load Bactopia/3.1.0

# Run Bactopia pipeline with TB-Profiler and spoligotyping enabled
bactopia \
    --path "$WORKDIR" \
    --datasets "$DBDIR" \
    --outdir "$OUTDIR" \
    --singularity \
    --tools tb-profiler \
    --arguments "--spoligotype"

# Confirm completion
echo "Spoligotyping analysis with Bactopia completed."
