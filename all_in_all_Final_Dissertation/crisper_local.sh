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



# Simulate Illumina reads from the genome
#wgsim -1 150 -2 150 -N 1000000 $GENOME $READS_DIR/P04_R1.fq $READS_DIR/P04_R2.fq
cd $HOME/CRISPRbuilder-TB
# Run CRISPRbuilder on the simulated reads
#python $HOME/CRISPRbuilder-TB/crisprbuilder.py \
#    -sra $READS_DIR \
#    -out $OUT_DIR \
#    -num_threads 12

python crisprbuilder.py \
  -sra /home/ma95362/crisprbuilder_test/P4_input \
  -out /home/ma95362/crisprbuilder_test/P4_output \
  -num_threads 12

