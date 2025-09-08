#!/bin/bash
#SBATCH --job-name=kraken2_job
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180gb
#SBATCH --time=05-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%A_%a.out
#SBATCH --error=/scratch/ma95362/scratch/log.%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu
#SBATCH --array=0-63   # adjust based on how many samples you have

#-----------------------------------
# Load Kraken2
#-----------------------------------
module load Kraken2/2.1.3-gompi-2023a

#-----------------------------------
# Input FASTQ directory
#-----------------------------------
FASTQ_DIR=/scratch/ma95362/ETH_bovis_Sequence/reads/64_SH_Sequence_data/raw/fixed_fastqs/fixed_fastqs
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
