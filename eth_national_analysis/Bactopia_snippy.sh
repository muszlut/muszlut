#!/bin/bash
#SBATCH --job-name=bactopia_snippy
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64          # use as many as your cluster allows
#SBATCH --mem=240G                  # bump memory if possible
#SBATCH --time=04-00:00:00          # 4 days walltime
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err            # Standard error log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# -------------------------
# Load module
# -------------------------
module load Bactopia/3.2.0

# -------------------------
# Paths
# -------------------------
RESULTS=/scratch/ma95362/eth_national_analysis/bactopia_results
SNIPPY_OUT=/scratch/ma95362/eth_national_analysis/snippy_results
REF="/scratch/ma95362/gbk/ncbi_dataset/data/GCF_000195955.2/genomic.gbk"

mkdir -p "$SNIPPY_OUT"
cd "$SNIPPY_OUT" || { echo "Failed to cd into $SNIPPY_OUT"; exit 1; }   
# -------------------------
# Run Snippy workflow
# -------------------------
bactopia --wf snippy \
  --bactopia "$RESULTS" \
  --outdir "$SNIPPY_OUT" \
  --cpus 64 \
  --reference "$REF"

echo "Snippy workflow finished!"
