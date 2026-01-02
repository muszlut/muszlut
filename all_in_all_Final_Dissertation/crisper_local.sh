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

cd /home/ma95362/CRISPRbuilder-TB

python crisprbuilder.py \
  -sra /home/ma95362/crisprbuilder_test/P4_readss/P04_combined.fastq.gz \
  -out /home/ma95362/crisprbuilder_test/P4_readss/output_P04 \
  -num_threads 12
