#!/bin/bash
#SBATCH --job-name=magma
#SBATCH --output=magma_%j.out
#SBATCH --error=magma_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=compute

# Load modules
module load Java/17.0.7 # or Java/11 if that’s what’s available
module load Miniconda3/23.5.2-0
module load Nextflow/23.10.0 # adjust to the available version
module load Mamba/23.1.0-4   # optional, if you use mamba

# (optional) Activate base conda if needed
# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate base

# Navigate to your working directory
cd /scratch/ma95362/magma_project

# Run the pipeline
nextflow run https://github.com/TORCH-Consortium/MAGMA \
  -profile conda_local,pbs \
  -r v1.1.1 \
  -params-file params.yml

# Notes:
# - adjust `-profile` if you use docker/podman instead of conda
# - adjust `-r v1.1.1` to the release version/tag you want
