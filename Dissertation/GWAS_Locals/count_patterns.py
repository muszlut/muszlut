#!/bin/bash
#SBATCH --job-name=count_patterns
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/ma95362/scratch/count_patterns_%j.out
#SBATCH --error=/scratch/ma95362/scratch/count_patterns_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# ------------------------------
# 1. Activate Conda environment
# ------------------------------
module load Miniforge3
source activate pyseer-env

# ------------------------------
# 2. Define paths
# ------------------------------
PYSEER_SCRIPT=/home/ma95362/pyseer/scripts/count_patterns.py
PATTERN_FILE=./pyseer_out/gene_patterns_AZM.txt
OUTPUT_FILE=./pyseer_out/gene_patterns_AZM_counts.txt

# ------------------------------
# 3. Run count_patterns.py
# ------------------------------
echo "Running count_patterns.py on $(basename $PATTERN_FILE)..."

python ${PYSEER_SCRIPT} ${PATTERN_FILE} > ${OUTPUT_FILE}

echo "✅ Pattern counting completed successfully."
echo "Output file saved to: ${OUTPUT_FILE}"

# ------------------------------
# 4. Deactivate Conda environment
# ------------------------------
conda deactivate
