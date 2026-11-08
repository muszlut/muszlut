#!/bin/bash
#SBATCH --job-name=01_setup_env
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# One-time setup: conda environment for the tree and GWAS steps
set -euo pipefail
module load Mamba/23.11.0-0

if conda env list | grep -q "^mtb_gwas "; then
    echo "mtb_gwas already exists - skipping creation"
else
    mamba create -y -n mtb_gwas -c conda-forge -c bioconda \
        python=3.11 snp-sites iqtree pyseer panaroo bcftools samtools mosdepth seqkit
fi

source activate mtb_gwas
echo "== versions (record these for the Methods)"
snp-sites -V 2>&1 | head -1 || true
(iqtree2 --version 2>&1 || iqtree --version 2>&1) | head -1 || true
pyseer --version 2>&1 || true
panaroo --version 2>&1 || true
python --version
 