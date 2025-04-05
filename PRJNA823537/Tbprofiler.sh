#!/bin/bash
#SBATCH --job-name=tb_profile_module_on_Sapelo2                # Job name
#SBATCH --partition=batch                                      # Partition (queue) name
#SBATCH --ntasks=1                                             # Run on a single CPU
#SBATCH --cpus-per-task=8                                      # Number of cores per task
#SBATCH --mem=40gb                                             # Job memory request
#SBATCH --time=07-00:00:00                                     # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err            # Standard error log

# Set input directory where fastq files are located
INPUT_DIR="/scratch/ma95362/PRJNA823537_ET125/Tbprofiler"
# Set output directory for TBProfiler results
OUTDIR="/scratch/ma95362/PRJNA823537_ET125/Tbprofiler_module"

# Create output directory if it doesn't exist
mkdir -p $OUTDIR

# Load TBProfiler module
module load TBProfiler/6.6.2

# Move to output directory
cd $OUTDIR

# Set custom directory for TBProfiler to use for databases
DB_DIR="/scratch/ma95362/PRJNA823537_ET125/TBProfiler_db"

# Create the directory if it doesn't exist
mkdir -p $DB_DIR

# Loop through your 17 isolates
for R1 in ${INPUT_DIR}/*_R1.fastq.gz; do
    base=$(basename "$R1" _R1.fastq.gz)
    R2="${INPUT_DIR}/${base}_R2.fastq.gz"

    echo "Running TBProfiler for: $base"

    tb-profiler profile \
        --read1 "$R1" \
        --read2 "$R2" \
        --prefix "$base" \
        --spoligotype \
        --threads 4
        --db_dir $DB_DIR

    echo "Finished: $base"
    echo "-------------------------------"
done
