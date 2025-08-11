#!/bin/bash
#SBATCH --job-name=All_Bactopia_TBprofiler                     # Job name
#SBATCH --partition=batch                                      # Partition (queue) name
#SBATCH --ntasks=1                                             # Run on a single CPU task
#SBATCH --cpus-per-task=8                                      # Number of CPU cores per task
#SBATCH --mem=40gb                                             # Memory per node
#SBATCH --time=07-00:00:00                                     # Time limit (7 days)
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # STDOUT log file
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err            # STDERR log file

#SBATCH --mail-type=END,FAIL                                   # Email notifications for job done & fail
#SBATCH --mail-user=ma95362@uga.edu                            # Email recipient

# ---------------- SETTINGS ----------------
OUTDIR="/scratch/ma95362/musse_MGA/all_reads_Bactopia_Analysis"

# Create output directory if it doesn't exist
if [ ! -d "$OUTDIR" ]; then
    echo "Creating output directory: $OUTDIR"
    mkdir -p "$OUTDIR"
fi

# Load required modules
module load Bactopia/3.2.0
module load TB-Profiler/6.6.5

# Move to the output directory
cd "$OUTDIR" || { echo "ERROR: Cannot change to directory $OUTDIR"; exit 1; }

# Run TB-Profiler on existing Bactopia results
echo "Starting TB-Profiler on Bactopia results in $OUTDIR at $(date)"

bactopia tools tb-profiler \
    --bactopia "$OUTDIR" \
    --threads "$SLURM_CPUS_PER_TASK" || { echo "ERROR: TB-Profiler failed"; exit 1; }

echo "TB-Profiler completed successfully at $(date)"
