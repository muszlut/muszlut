#!/bin/bash
#SBATCH --job-name=download_H37Rv
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu  

OUTDIR="/scratch/ma95362/ggcaller_module/reference"

if [ ! -d $OUTDIR ]; then 
    mkdir -p $OUTDIR
fi

cd $OUTDIR

# Correct URL (no spaces!)
URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz"

# Download H37Rv
wget $URL

# Unzip
gunzip GCF_000195955.2_ASM19595v2_genomic.fna.gz
