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
#SBATCH --array=1-6

# Load environment
source ~/.bashrc
conda activate mtbvartools

# Paths
OUTPUT_DIR="/scratch/ma95362/Sequence"
SCRIPT_PATH="/scratch/ma95362/mtbvartools/scripts/sra_download.py"  # this is the script path to download SRA sequences from mtbvartools/scripts

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

SRA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sra_list.txt)
OUTPUT="sample_${SRA}"

python "$SCRIPT_PATH" \
    -i $SRA \
    -o $OUTPUT \
    -d ./downloads \
    --tmp-path ./tmp \
    --overwrite
