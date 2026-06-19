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

export TMPDIR=/scratch/ma95362/tmp
mkdir -p $TMPDIR

singularity exec /apps/singularity-images/ggcallaroo_v0.1.0.sif \
snakemake \
-s /app/Snakefile \
--directory /scratch/ma95362/257_assembled_files/MTB_ggcaller/test-ggcallaroo_results \
--cores 32 \
--config \
refs=/scratch/ma95362/257_assembled_files/MTB_ggcaller/test-ggcallaroo_results/test_refs.txt \
bakta_db=/scratch/ma95362/full_bakta_db/db \
output_dir=/scratch/ma95362/257_assembled_files/MTB_ggcaller/test-ggcallaroo_results \
ggcaller_cli_args="--save --balrog-db /scratch/ma95362/ggcaller_db" \
panaroo_cli_args="--clean-mode moderate"