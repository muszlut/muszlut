#!/bin/bash
#SBATCH --job-name=TBprofiler_Run                    # Job name
#SBATCH --partition=batch                            # Queue/partition
#SBATCH --ntasks=1                                  # Number of tasks (1 for single-process)
#SBATCH --cpus-per-task=8                           # Number of CPU cores
#SBATCH --mem=40gb                                  # Memory allocation
#SBATCH --time=07-00:00:00                          # Max runtime (7 days)
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out  # STDOUT log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err   # STDERR log

#SBATCH --mail-type=END,FAIL                         # Email notifications
#SBATCH --mail-user=ma95362@uga.edu                  # Email recipient

# Load required modules
module load Bactopia/3.2.0
module load TB-Profiler/6.6.5

# Directory where sample file is located
SAMPLE_DIR="/scratch/ma95362/musse_MGA/all_reads_Bactopia_Analysis"
SAMPLES_FILE="${SAMPLE_DIR}/samples_clean.fofn"

# Check if samples.fofn exists
if [ ! -f "$SAMPLES_FILE" ]; then
    echo "ERROR: Samples file $SAMPLES_FILE not found!"
    exit 1
fi

cd "$SAMPLE_DIR"

echo "Running TB-Profiler on samples listed in $SAMPLES_FILE"
echo "Job started at: $(date)"

# Run TB-Profiler with specified samples and threads
bactopia tools tb-profiler --samples "$SAMPLES_FILE" --threads "$SLURM_CPUS_PER_TASK" || { echo "ERROR: TB-Profiler failed"; exit 1; }

echo "TB-Profiler completed successfully at: $(date)"
