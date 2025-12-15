#!/bin/bash
#SBATCH --job-name=gunzip_job          # Job name
#SBATCH --partition=batch
#SBATCH --ntasks=1                     # One task
#SBATCH --cpus-per-task=4              # Number of CPU cores
#SBATCH --mem=8G                       # Memory (adjust as needed)
#SBATCH --time=02:00:00                # Max runtime (hh:mm:ss)
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# -----------------------------
# Load any required modules
# -----------------------------
#module load gzip/1.12   # or module load gzip (if available)
#cd /scratch/ma95362/all_in_all_reads/bactopia_prepare/bactopia-runs/pangenome-20251128-070449
cd /scratch/ma95362/new_ggcaller_pro
# -----------------------------
# Define input and output files
# -----------------------------
# Example 1: decompress a single file
gunzip -c ERR6397190.ffn.gz > ERR6397190.ffn
gunzip -c ERR6397187.ffn.gz > ERR6397187.ffn
gunzip -c ERR6397199.ffn.gz > ERR6397199.ffn
gunzip -c ERR6397200.ffn.gz > ERR6397200.ffn