#!/bin/bash
#SBATCH --job-name=ggcaller_new
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=40G
#SBATCH --time=10:00:00
#SBATCH --output=/scratch/ma95362/test_ggcaller_out/log.%j.out
#SBATCH --error=/scratch/ma95362/test_ggcaller_out/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Direct path to ggcaller binary inside your conda environment
GG=/home/ma95362/.conda/envs/ggc_env/bin/ggcaller

OUTDIR=test_ggcaller_out
mkdir -p $OUTDIR

$GG \
    --reads1 reads_1.fq.gz \
    --reads2 reads_2.fq.gz \
    --outdir $OUTDIR \
    --threads 4
