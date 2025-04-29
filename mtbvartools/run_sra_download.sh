#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH --partition=batch                             # Partition (queue) name
#SBATCH --ntasks=1                                    # Run on a single CPU
#SBATCH --cpus-per-task=8                             # Number of cores per task
#SBATCH --mem=40GB                                    # Job memory request
#SBATCH --time=07-00:00:00                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out  # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err   # Standard error log

#SBATCH --mail-type=END,FAIL                          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                   # Where to send mail

module purge
module load python/3.10

# (Optional) activate your environment
source activate mtbvartools

# Create logs directory if not exists
mkdir -p logs

# Read the nth line of the SRA list (based on SLURM_ARRAY_TASK_ID)
SRA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sra_list.txt)

# Set output name using the accession
OUTPUT="sample_${SRA}"

# Run your Python script
python sra_download_script.py \
    -i $SRA \
    -o $OUTPUT \
    -d ./downloads \
    --tmp-path ./tmp \
    --overwrite
