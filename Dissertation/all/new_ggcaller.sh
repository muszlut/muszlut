#!/bin/bash
#SBATCH --job-name=ggcaller_new
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=40G
#SBATCH --time=10:00:00
#SBATCH --output=/scratch/ma95362/test_ggcaller_out/log.%j.out
#SBATCH --error=/scratch/ma95362/test_ggcaller_out/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

cd /scratch/ma95362/ggcaller_reads

# Load modules or environment
module purge
module load miniconda3          # Or appropriate module on your cluster

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ggc_env
#bifrost build \
#    -r test_reads.txt \
#    -k 31 \
#    -o bifrost_graph_reads \
#    -t 8 \
#    --colors



# Path to graph and output
ggcaller \
    --graph bifrost_graph_reads.gfa \
    --colours bifrost_graph_reads.bfg_colors \
    --out ggcaller_output \
    --threads 8 \
    --kmer 31

# Optional: print completion message
echo "ggCaller run finished at $(date)"
