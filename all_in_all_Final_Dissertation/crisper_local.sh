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

cd /home/ma95362/crisprbuilder_test/P4_readss

python /home/ma95362/CRISPRbuilder-TB/crisprbuilder.py \
  -i P04._R1.fastq.gz \
  -j P04._R2.fastq.gz \
  -o output_P04
