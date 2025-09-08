#!/bin/bash
#SBATCH --job-name=kraken2_summary
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Load modules
module load Kraken2/2.1.3-gompi-2023a
module load Python/3.11.5

# Path to your KrakenTools installation
KRAKENTOOLS=/scratch/ma95362/KrakenTools   # <-- update if needed
REPORTS_DIR=/scratch/ma95362/ETH_bovis_Sequence/reads/64_SH_Sequence_data/raw/fixed_fastqs/fixed_fastqs
OUTFILE=/scratch/ma95362/ETH_bovis_Sequence/reads/64_SH_Sequence_data/raw/fixed_fastqs/fixed_fastqs/kraken2_combined.tsv

# Combine reports into a single abundance table
python $KRAKENTOOLS/combine_kreports.py \
    -r ${REPORTS_DIR}/*.txt \
    -o ${OUTFILE} \
    -p
