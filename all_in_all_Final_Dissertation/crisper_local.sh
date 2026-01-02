#!/bin/bash
#SBATCH --job-name=crispr_P04
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=36:00:00
#SBATCH --output=/scratch/ma95362/crisprbuilder_test/logs/log.%j.out
#SBATCH --error=/scratch/ma95362/crisprbuilder_test/logs/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Load conda
source ~/.bashrc
conda activate crisprbuilder_tb

# Define paths
GENOME=/home/ma95362/crisprbuilder_test/P4/P04.fna
READS_DIR=/home/ma95362/crisprbuilder_test/P4_reads
OUT_DIR=/home/ma95362/crisprbuilder_test/P4_out

# Create reads directory
mkdir -p $READS_DIR
mkdir -p $OUT_DIR

# Simulate Illumina reads from the genome
wgsim -1 150 -2 150 -N 1000000 $GENOME $READS_DIR/P04_R1.fq $READS_DIR/P04_R2.fq

# Run CRISPRbuilder on the simulated reads
python $HOME/CRISPRbuilder-TB/crisprbuilder.py \
    -sra $READS_DIR \
    -out $OUT_DIR \
    -num_threads 12
