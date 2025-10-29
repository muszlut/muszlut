#!/bin/bash
#SBATCH --job-name=R_eptb_ptb
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/ma95362/scratch/R_eptb_ptb_%j.out
#SBATCH --error=/scratch/ma95362/scratch/R_eptb_ptb_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

#Set output directory variable
OUTDIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_local_all/pyseer_output"
#Location of R script
TIDY="/home/ma95362/muszlut/Dissertation/GWAS_Locals/pyseer_anlaysis_EPTB_PTB.R"


#Tell the program to make  the outdir folder
if [ ! -d $OUTDIR ] 
    then 
        mkdir -p $OUTDIR
fi

cd $OUTDIR

#Activate Renvironmemt previously created
source activate r-tidyverse

#Run R script:here
R --no-save < $TIDY
