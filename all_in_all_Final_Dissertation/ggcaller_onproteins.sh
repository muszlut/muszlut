#!/bin/bash
#SBATCH --job-name=ggcaller_somali_only
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=895G
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/ma95362/ggcaller_module/logs/somali/log.%j.out
#SBATCH --error=/scratch/ma95362/ggcaller_module/logs/somali/log.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

module load ggCaller/1.4.1
cd /scratch/ma95362/new_ggcaller_pro

# scalable input
ls ggcaller_input/*.ffn > input.txt

ggcaller \
  --refs input.txt \
  --identity-cutoff 0.95 \
  --len-diff-cutoff 0.9 \
  --family-threshold 0.7 \
  --alignment pan \
  --aligner def \
  --out ggc_mtbc \
  --threads 8 \
    --balrog-db /scratch/ma95362/ggcaller_db/ggCallerdb