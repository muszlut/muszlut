#!/bin/bash
#SBATCH --job-name=CRISPRbuilder_SRA
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/ma95362/CRISPRbuilder-TB/logs.%A_%a.out
#SBATCH --error=/scratch/ma95362/CRISPRbuilder-TB/logs.%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu
#SBATCH --array=1-257%20 

# -------------------------
# 1. Environment & Paths
# -------------------------
module purge
source ~/.bashrc
conda activate crisprbuilder_tb

BASE_DIR=/scratch/ma95362/CRISPRbuilder-TB
SAMPLE_FILE=$BASE_DIR/samples.txt  # The list of sample names

# Move to the project directory
cd $BASE_DIR || exit 1

# -------------------------
# 2. Get Sample Name
# -------------------------
# Extract the sample name from the specific line in samples.txt matching the Array ID
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLE_FILE)

echo "Processing Task ID: $SLURM_ARRAY_TASK_ID"
echo "Sample Name: $SAMPLE_NAME"
echo "Working Directory: $(pwd)"

# -------------------------
# 3. Run CRISPRbuilder
# -------------------------
# Use the sample name (folder name) as the -sra argument.
# The tool will look inside the 'sequences/' directory for this name.
python crisprbuilder.py -sra "$SAMPLE_NAME" -num_threads 12