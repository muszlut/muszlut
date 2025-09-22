#!/bin/bash
#SBATCH --job-name=bovisanalyzer
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=72:00:00
#SBATCH --mem=128G
#SBATCH --output=bovisanalyzer_%j.out
#SBATCH --error=bovisanalyzer_%j.err

# Load Anaconda / Miniconda
module load anaconda/2020.11

# Activate the bovisanalyzer environment
source activate bovisanalyzer

# Define paths
SAMPLESHEET=/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples/samplesheet.csv
REFERENCE=/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/AF2122_97.fasta
KRAKEN2DB=/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/minikraken2_v1_8GB
OUTDIR=/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_results

# Make output directory if it doesn't exist
mkdir -p $OUTDIR

# Run Bovisanalyzer
nextflow run avantonder/bovisanalyzer \
    -profile conda \
    --input $SAMPLESHEET \
    --reference $REFERENCE \
    --kraken2db $KRAKEN2DB \
    --outdir $OUTDIR
