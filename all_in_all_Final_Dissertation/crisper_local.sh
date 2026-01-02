#!/bin/bash
#SBATCH --job-name=SPAdes_CRISPR
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/ma95362/crisprbuilder_test/logs/log.%j.out
#SBATCH --error=/scratch/ma95362/crisprbuilder_test/logs/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# ===============================
# Activate conda environment
# ===============================
source ~/miniconda3/etc/profile.d/conda.sh
conda activate crisprbuilder_tb

# ===============================
# Fixed sample paths (YOUR DATA)
# ===============================
SAMPLE="P01"
ASSEMBLY="/home/ma95362/crisprbuilder_test/P01/P01_spades/contigs.fasta"
CRISPR_DIR="$HOME/CRISPRbuilder-TB"

# ===============================
# Prepare directories
# ===============================
mkdir -p "$CRISPR_DIR/data"
mkdir -p "$CRISPR_DIR/results"
mkdir -p logs

# ===============================
# Copy contigs into CRISPRbuilder
# ===============================
cp "$ASSEMBLY" "$CRISPR_DIR/data/${SAMPLE}.fasta"

# ===============================
# Run CRISPRbuilder-TB
# ===============================
cd "$CRISPR_DIR" || exit 1

python CRISPRbuilder.py \
  -i "data/${SAMPLE}.fasta" \
  -o "results/${SAMPLE}"

echo "✔ CRISPRbuilder completed successfully for ${SAMPLE}"