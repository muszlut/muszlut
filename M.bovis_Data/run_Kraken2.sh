#!/bin/bash
#SBATCH --job-name=Bactopia_Kraken2                     # Job name
#SBATCH --partition=batch                                # Partition (queue)
#SBATCH --ntasks=1                                       # Single task
#SBATCH --cpus-per-task=8                                # CPUs per task
#SBATCH --mem=40gb                                       # Memory
#SBATCH --time=05-00:00:00                               # Time limit (HH:MM:SS)
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out    # STDOUT log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err     # STDERR log
#SBATCH --mail-type=END,FAIL                             # Mail events
#SBATCH --mail-user=ma95362@uga.edu                      # Your email


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

