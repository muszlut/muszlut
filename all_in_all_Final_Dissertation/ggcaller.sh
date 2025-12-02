#!/bin/bash
#SBATCH --job-name=ggcaller_CPS
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/ma95362/ggcaller_module/logs/log.%j.out
#SBATCH --error=/scratch/ma95362/ggcaller_module/logs/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

echo "Job started on $(date)"

# Load the ggcaller module installed by GACRC
module load ggCaller/1.4.1

# Move to the directory where you submitted the job
cd /scratch/ma95362/ggcaller_module

# Run ggCaller
ggcaller --refs Bentley_et_al_2006_CPS_sequences/input.txt \
         --annotation ultrasensitive \
         --diamonddb Bentley_et_al_2006_CPS_protein_sequences.faa \
         --aligner def \
         --alignment pan \
         --save \
         --out ggc_Bentley_et_al_CPS \
         --threads 4 \
         --balrog-db /scratch/ma95362/

echo "Job finished on $(date)"
