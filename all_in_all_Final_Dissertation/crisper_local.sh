#!/bin/bash
#SBATCH --job-name=SPAdes_CRISPR
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/ma95362/crisprbuilder_test/logs/log.%j.out
#SBATCH --error=/scratch/ma95362/crisprbuilder_test/logs/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# ===============================
# Activate conda environment
# ===============================
source ~/.bashrc
conda activate crisprbuilder_tb

# Directories
CRISPR_DIR="$HOME/CRISPRbuilder-TB"
mkdir -p "$CRISPR_DIR/data" "$CRISPR_DIR/results" logs

# Make sure your FASTA is in the data folder
# cp /path/to/P01/contigs.fasta $CRISPR_DIR/data/P01.fasta

cd "$CRISPR_DIR" || exit 1

# Run CRISPRbuilder-TB on all FASTA files in data/
python crisprbuilder.py -out results/P01 -num_threads 12

echo "✔ CRISPRbuilder completed successfully for P01"
