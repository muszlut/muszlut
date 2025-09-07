#!/bin/bash
#SBATCH --job-name=kraken2_bactopia
#SBATCH --partition=batch
#SBATCH --output=kraken2_%A_%a.log
#SBATCH --error=kraken2_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=05-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu
#SBATCH --array=0-63   # adjust to total number of samples minus 1

#-----------------------------------
# Load Bactopia (Kraken2 included)
#-----------------------------------
module load Bactopia/3.2.0

#-----------------------------------
# Input FASTQ directory
#-----------------------------------
FASTQ_DIR=/scratch/ma95362/ETH_bovis_Sequence/reads/64_SH_Sequence_data/raw/fixed_fastqs
cd $FASTQ_DIR || { echo "Directory not found: $FASTQ_DIR"; exit 1; }

#-----------------------------------
# Sample list
#-----------------------------------
SAMPLES=($(ls *_R1.fastq.gz | sort))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
BASE=$(basename $SAMPLE _R1.fastq.gz)

#-----------------------------------
# Kraken2 database
#-----------------------------------
KRAKEN2_DB=/scratch/ma95362/kraken2_db

#-----------------------------------
# Run Kraken2
#-----------------------------------
kraken2 \
  --db $KRAKEN2_DB \
  --paired ${BASE}_R1.fastq.gz ${BASE}_R2.fastq.gz \
  --threads 8 \
  --report ${BASE}_kraken2_report.txt \
  --output ${BASE}_kraken2_output.txt

