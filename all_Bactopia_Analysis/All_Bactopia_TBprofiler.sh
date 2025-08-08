#!/bin/bash
#SBATCH --job-name=All_Bactopia_TBprofiler                     # Job name
#SBATCH --partition=batch                                      # Partition (queue) name
#SBATCH --ntasks=1                                             # Run on a single CPU
#SBATCH --cpus-per-task=8                                      # Number of cores per task
#SBATCH --mem=40gb                                             # Job memory request
#SBATCH --time=07-00:00:00                                     # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err            # Standard error log

#SBATCH --mail-type=END,FAIL                                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                            # Where to send mail	

# ---- SETTINGS ----
BACTOPIA_OUT="/scratch/ma95362/musse_MGA/all_reads_Bactopia_Analysis/MGA_paired_end_samples"
TBPROFILER_OUT="/scratch/ma95362/musse_MGA/all_reads_Bactopia_Analysis/MGA_tbprofiler_results"

# Create output directory if needed
mkdir -p $TBPROFILER_OUT

# Load modules
module load Bactopia/3.2.0
module load TB-Profiler/6.6.5

cd $TBPROFILER_OUT
# ---- STEP 1: Run TB-Profiler via Bactopia tools ----
echo "Running TB-Profiler on existing Bactopia results..."
bactopia-tools tb-profiler \
    --bactopia $BACTOPIA_OUT \
    --include "E*.,P*." \
    --outdir $TBPROFILER_OUT \
    --threads $SLURM_CPUS_PER_TASK


# ---- STEP 2: Collate TB-Profiler results ----
echo "Collating TB-Profiler results..."
cd $TBPROFILER_OUT
tb-profiler collate
cd -

echo "All done."
echo "TB-Profiler results are in: $TBPROFILER_OUT"