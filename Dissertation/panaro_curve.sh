#!/bin/bash
#SBATCH --job-name=panaroo_curve
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
OUTPUT_DIR="$INPUT_DIR/Plot_output"


# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Move to the input directory
cd "$INPUT_DIR"

# Run Panaroo filter
panaroo-plot --input /scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_of_1368/panaroo \
             --plot-type rarefaction \
             --output panaroo_rarefaction_curve.png

echo "✅ Panaroo plot completed successfully."

conda deactivate