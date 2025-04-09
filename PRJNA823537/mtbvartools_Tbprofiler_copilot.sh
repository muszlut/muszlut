#!/bin/bash
#SBATCH --job-name=mtbvarTools_tbprofiler_copilot_job # Job name
#SBATCH --partition=batch                             # Partition (queue) name
#SBATCH --ntasks=1                                    # Run on a single CPU
#SBATCH --cpus-per-task=8                             # Number of cores per task
#SBATCH --mem=40gb                                    # Job memory request
#SBATCH --time=07-00:00:00                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out  # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err   # Standard error log

#SBATCH --mail-type=END,FAIL                          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                   # Where to send mail

# Activate the Conda environment
source ~/.bashrc
conda activate mtbvartools

# Define directories
BAM_DIR="/scratch/ma95362/PRJNA823537_ET125/mtbVARTolls_trial"
OUTPUT_DIR="/scratch/ma95362/PRJNA823537_ET125/mtbVARTolls_trial/tbprofiler_copilot"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Change to output directory
cd "$OUTPUT_DIR"

# Load Python in SLURM session
python <<EOF
import os
import mtbvartools as vt

# Get all BAM files in the directory
bam_files = [f for f in os.listdir("$BAM_DIR") if f.endswith(".bam")]

for bam in bam_files:
    sample_name = bam.replace(".bam", "")
    bam_path = os.path.join("$BAM_DIR", bam)
    output_prefix = os.path.join("$OUTPUT_DIR", sample_name)

    print(f"Processing sample: {sample_name}")

    # Run TBProfiler through MTBVarTools
    cmd = f"tb-profiler profile --bam {bam_path} --threads 8 --prefix {output_prefix}"
    vt.contShell(cmd)

print("TBProfiler analysis complete!")
EOF

echo "Job finished!"