#!/bin/bash
#SBATCH --job-name=ggcaller_graph
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/ma95362/test_ggcaller_out/log.%j.out
#SBATCH --error=/scratch/ma95362/test_ggcaller_out/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

echo "START: $(date)"

# -------------------------
# Load mamba on Sapelo2
# -------------------------
module load Mamba/23.11.0-0

# Activate your environment
source activate ~/.conda/envs/ggc_env
# OR (if above fails):
conda activate ggc_env

# -------------------------
# Move to submission folder
# -------------------------
cd /scratch/ma95362/ggcaller_reads
# -------------------------
# 1. Build colored Bifrost graph
# -------------------------
#bifrost build \
#  -r test_reads.txt \
#  -k 31 \
#  -o bifrost_graph_reads \
#  --colors \
#  --threads 16

# -------------------------
# 2. Run ggCaller with color file
# -------------------------
ggcaller \
  --graph test_reads.gfa \
  --colours test_reads.bfg_colors \
  --reads-list test_reads.txt \
  --kmer 31 \
  --threads 16 \
  --out ggCaller_output

echo "END: $(date)"
