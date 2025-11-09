#!/bin/bash
#SBATCH --job-name=Bactopia_Full_Run              # Job name
#SBATCH --partition=batch                         # Partition name
#SBATCH --ntasks=1                                # Number of tasks
#SBATCH --cpus-per-task=16                        # Number of CPU cores
#SBATCH --mem=180gb                               # Total memory
#SBATCH --time=03-00:00:00                        # Time limit (3 days)
#SBATCH --output=/scratch/ma95362/logs/bactopia_full_run.%j.out   # Standard output
#SBATCH --error=/scratch/ma95362/logs/bactopia_full_run.%j.err    # Standard error
#SBATCH --mail-type=END,FAIL                      # Mail when job ends or fails
#SBATCH --mail-user=ma95362@uga.edu               # Your email

# -------------------------------
# 1️⃣  Load the environment
# -------------------------------
module load Bactopia/3.2.0-conda

# -------------------------------
# 2️⃣  Define paths
# -------------------------------
WORKDIR="/scratch/ma95362/EPTB_Hilina/ETH_Bactopia_Prepare"
OUTDIR="${WORKDIR}/ETH_full_analysis"
SAMPLES="${WORKDIR}/ETH_samples.txt"

# Make sure output directories exist
mkdir -p "$OUTDIR"
mkdir -p /scratch/ma95362/logs

# -------------------------------
# 3️⃣  Run the full Bactopia pipeline
# -------------------------------
bactopia \
  --samples "$SAMPLES" \
  --species "Mycobacterium tuberculosis" \
  --genome-size 4410000 \
  --coverage 100 \
  --datasets /apps/eb/Bactopia/3.2.0/share/bactopia-3.2.0/datasets/ \
  --outdir "$OUTDIR" \
  --max_cpus 16 \
  --max_time '72.h'

# -------------------------------
# 4️⃣  Generate a summary report
# -------------------------------
bactopia summary \
  --bactopia-path "$OUTDIR"
