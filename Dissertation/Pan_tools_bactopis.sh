#!/bin/bash
#SBATCH --job-name=BacPanTool_new
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# -----------------------------
# Load environment modules
# -----------------------------
module load Bactopia/3.2.0-conda
module load Java/17.0.6

# -----------------------------
# Set directories
# -----------------------------
BACTOPIA_DIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples"
OUTDIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/bactopia_Pangenome_tools_results"

# -----------------------------
# Create output folder if not exist
# -----------------------------
mkdir -p $OUTDIR
cd $OUTDIR

# -----------------------------
# Run Bactopia pangenome tool
# -----------------------------
bactopia tools pangenome \
    --dir $BACTOPIA_DIR \
    --exclude $BACTOPIA_DIR/bactopia-exclude_final.tsv \
    --outdir $OUTDIR
