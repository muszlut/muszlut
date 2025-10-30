#!/bin/bash
#SBATCH --job-name=panaroo_filter_pa
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=07-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Load conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate panaroo-env

# Define directories
INPUT_DIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_of_1368/panaroo"
INPUT_FILE="gene_presence_absence.csv"  # just the filename!
OUTPUT_DIR="$INPUT_DIR/filtered_output"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Move to the input directory
cd "$INPUT_DIR"

# Run Panaroo filter
panaroo-filter-pa --input "$INPUT_FILE" --out_dir "$OUTPUT_DIR" --type pseudo,length

echo "✅ Panaroo filter-pa completed successfully."

conda deactivate