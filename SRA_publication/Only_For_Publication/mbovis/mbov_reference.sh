#!/bin/bash
#SBATCH --job-name=get_mbovis_ref
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --mem=2GB
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/ma95362/ref_mbovis/log.%j.out
#SBATCH --error=/scratch/ma95362/ref_mbovis/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

set -euo pipefail

OUTDIR="/scratch/ma95362/ref_mbovis"
ACCESSION="GCF_000195835.3"     # M. bovis AF2122/97

mkdir -p "$OUTDIR"
cd "$OUTDIR"

module load NCBI-Datasets-CLI/16.4.4

datasets download genome accession "$ACCESSION" \
    --include genome,gbff,cds,protein,gff3 \
    --filename mbovis_af2122.zip

unzip -o mbovis_af2122.zip

# show what arrived
find ncbi_dataset -type f | sed 's/^/  /'