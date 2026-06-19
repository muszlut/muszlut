#!/bin/bash
#SBATCH --job-name=ggCallaroo
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120gb
#SBATCH --time=07-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

snakemake \
-s /app/Snakefile \
--directory /scratch/ma95362/257_assembled_files/MTB_ggcaller/test-ggcallaroo_results \
--cores 32 \
--use-conda \
--dry-run \
--config \
refs=/scratch/ma95362/257_assembled_files/MTB_ggcaller/test-ggcallaroo_results/test_refs.txt \
bakta_db=/scratch/ma95362/full_bakta_db/db \
output_dir=/scratch/ma95362/257_assembled_files/MTB_ggcaller/test-ggcallaroo_results \
ggcaller_cli_args="--save --balrog-db /scratch/ma95362/ggcaller_db" \
panaroo_cli_args="--clean-mode moderate"