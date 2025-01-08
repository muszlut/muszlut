#!/bin/bash
#SBATCH --job-name=Musse_WGS_collate           # Job name
#SBATCH --partition=batch                      # Partition (queue) name
#SBATCH --ntasks=1                             # Run on a single CPU
#SBATCH --cpus-per-task=2                      # Number of cores per task
#SBATCH --mem=10gb                             # Job memory request
#SBATCH --time=01:00:00                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out   # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu            # Where to send mail

# Set output directory variable
OUTDIR="/scratch/ma95362/musse_MGA/fastqs/MGA_paired_end_samples"
# Define path to the sequence reads
READ1="/work/fdqlab/Ethiopia_wgs_mtb_2024/first_run/sample_read1.fastq"
READ2="/work/fdqlab/Ethiopia_wgs_mtb_2024/first_run/sample_read2.fastq"

# Create output directory if it doesn't exist
if [ ! -d "$OUTDIR" ]; then
    mkdir -p "$OUTDIR"
fi

# Load necessary modules
module load Miniconda3/23.5.2-0

# Activate conda environment (should be version 6.3.0)
source activate tb-profiler-env

# Move to working directory
cd "$OUTDIR"

# Run tb-profiler spoligotype analysis with paired-end reads
tb-profiler spoligotype --read1 "$READ1" --read2 "$READ2"

# Print a message indicating completion
echo "Spoligotyping analysis complete. Results are in $OUTDIR"
