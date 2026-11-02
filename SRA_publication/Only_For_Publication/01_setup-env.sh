#!/bin/bash
#SBATCH --job-name=01_setup_env
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --time=07-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# One-time setup: a conda environment with everything needed for the tree and GWAS steps.
# Run on an interactive node (not the login node):  interact --mem=8G  then  bash 01_setup_env.sh
set -euo pipefail
module load Mamba/23.11.0-0
mamba create -y -n mtb_gwas -c conda-forge -c bioconda \
    python=3.11 snp-sites iqtree pyseer panaroo bcftools samtools mosdepth seqkit
source activate mtb_gwas
echo "== versions (record these for the Methods)"
snp-sites -V 2>&1 | head -1
iqtree2 --version 2>&1 | head -1 || iqtree --version 2>&1 | head -1
pyseer --version 2>&1
panaroo --version 2>&1
 