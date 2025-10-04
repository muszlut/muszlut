#!/bin/bash
#SBATCH --job-name=bovis_analyzer                 # Job name
#SBATCH --partition=batch                         # Partition (queue)
#SBATCH --ntasks=1                                # Single task
#SBATCH --cpus-per-task=8                         # CPUs per task
#SBATCH --mem=180G                                # Memory
#SBATCH --time=05-00:00:00                        # Time limit (HH:MM:SS)
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out   # STDOUT log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err    # STDERR log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Initialize Conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bovisanalyzer

# Define paths
SAMPLESHEET=/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples/samplesheet.csv
REFERENCE=/scratch/ma95362/ETH_bovis_Sequence/bovis_REF/Fasta/AF2122_97.fasta
kraken2_db=/scratch/ma95362/kraken2_db
OUTDIR=/scratch/ma95362/ETH_bovis_Sequence/bovisanalyzer_output

# Make directories
mkdir -p $OUTDIR
mkdir -p $kraken2_db
cd $OUTDIR

# Run Bovisanalyzer
nextflow run avantonder/bovisanalyzer \
    -profile conda \
    --input $SAMPLESHEET \
    --reference $REFERENCE \
    --kraken2db $kraken2_db \
    --outdir $OUTDIR \
    -resume \
    -c /scratch/ma95362/bovisanalyzer_custom.config

