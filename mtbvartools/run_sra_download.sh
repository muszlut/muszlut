#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40GB
#SBATCH --time=07-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/logs/sra_download_%A_%a.out
#SBATCH --error=/scratch/ma95362/scratch/logs/sra_download_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu
#SBATCH --array=1-5

module purge
module load python/3.10

# Activate your environment
source activate mtbvartools

# Paths
OUTPUT_DIR="/scratch/ma95362/Sequence"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# Read current SRA ID
SRA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sra_list.txt)
OUTPUT="sample_${SRA}"

# Run script
python sra_download_script.py \
    -i $SRA \
    -o $OUTPUT \
    -d ./downloads \
    --tmp-path ./tmp \
    --overwrite
