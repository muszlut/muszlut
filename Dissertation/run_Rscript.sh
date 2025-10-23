#!/bin/bash
#SBATCH --job-name=Rscript_T3_ETH
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40gb
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu
#run
OUTDIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/pangenome_tools_results/bactopia/bactopia-runs/pangenome_of_L4.2.2.2.2/panaroo/filtered_output/pyseer_drugtype_out"
TIDY="/home/ma95362/muszlut/Dissertation/dr_gwas.R"

module load Miniforge3
source ~/.bashrc
conda activate r-tidyverse

cd $OUTDIR
Rscript $TIDY
