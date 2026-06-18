#!/bin/bash
#SBATCH --job-name=bakta_db
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/ma95362/bakta_logs/bakta_db_%j.out
#SBATCH --error=/scratch/ma95362/bakta_logs/bakta_db_%j.err

mkdir -p /scratch/ma95362/bakta_logs
mkdir -p /scratch/ma95362/bakta_db
cd /scratch/ma95362/bakta_db
singularity exec /apps/singularity-images/ggcallaroo_v0.1.0.sif \
bakta_db download \
    --output /scratch/ma95362/bakta_db \
    --type light