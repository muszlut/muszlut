#!/bin/bash
#SBATCH --job-name=crispr_SRR26800480
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

# Move to software directory (important for data/SIT.xls)
cd $HOME/CRISPRbuilder-TB

# Run CRISPRbuilder
#python crisprbuilder.py \
#  -sra /home/ma95362/crisprbuilder_test/SRR26800480 \
#  -out /home/ma95362/crisprbuilder_test/SRR26800480 \
#  -num_threads 12

python crisprbuilder.py \
  -sra /home/ma95362/crisprbuilder_test/P4 \
  -out /home/ma95362/crisprbuilder_test/P4 \
  -num_threads 12