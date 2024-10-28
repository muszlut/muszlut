#!/bin/bash
#SBATCH --job-name=Bactopia_tutorial        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=4              # Number of cores per task
#SBATCH --mem=40gb                     # Job memory request
#SBATCH --time=03:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err             # Standard error log

#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu  # Where to send mail	

#Set output directory variable
OUTDIR="/scratch/ma95362/PRJNA1056148_bactopia/ena-multiple-samples/bactopia-runs/snippy-20241021-142944/gubbins"
#Location of R script
TIDY="/home/ma95362/muszlut/coun_snps.R"


#Tell the program to make  the outdir folder
if [ ! -d $OUTDIR ] 
    then 
        mkdir -p $OUTDIR
fi

cd $OUTDIR

#gunzip summary of snp distribution file
gunzip -c +4 $OUTDIR/core-snp.summary_of_snp_distribution.vcf.gz > $OUTDIR/core-snp.summary_of_snp_distribution.vcf

#Remove header:
tail -n +4 $OUTDIR/core-snp.summary_of_snp_distribution.vcf > $OUTDIR/1core-snp.summary_of_snp_distribution.vcf

#Covert to csv
awk 'BEGIN {OFS=","} {$1=$1} 1' $OUTDIR/1core-snp.summary_of_snp_distribution.vcf > $OUTDIR/Ethiopia_core-snp.summary_of_snp_distribution.csv

#Activate Renvironmemt previously created
source activate r-tidyverse

#Run R script:here
R --no-save < $TIDY
