#!/bin/bash
#SBATCH --job-name=bovis_analyzer
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=05-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Load conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bovisanalyzer

# Define inputs
SAMPLESHEET=/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples/samplesheet.csv
REFERENCE=/scratch/ma95362/ETH_bovis_Sequence/bovis_REF/Fasta/AF2122_97.fasta
KRAKEN2DB=/scratch/ma95362/kraken2_db/mini_db
OUTDIR=/scratch/ma95362/ETH_bovis_Sequence/bovisanalyzer_output
#CONFIG=/scratch/ma95362/bovisanalyzer_custom.config
CONFIG=/scratch/ma95362/bovisanalyzer_resume.config

# Create output and logs directory
mkdir -p $OUTDIR
cd $OUTDIR

# Run Bovisanalyzer pipeline bovisanalyzer
nextflow run avantonder/bovisanalyzer \
    -c $CONFIG \
    --input $SAMPLESHEET \
    --reference $REFERENCE \
    --kraken2db $KRAKEN2DB \
    --outdir $OUTDIR \
    -resume
