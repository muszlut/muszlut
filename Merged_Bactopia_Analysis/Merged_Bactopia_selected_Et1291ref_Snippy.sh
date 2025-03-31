#!/bin/bash
#SBATCH --job-name=Musse_Merged_Et129ref_snippy_Excluded       # Job name
#SBATCH --partition=batch                                      # Partition (queue) name
#SBATCH --ntasks=1                                             # Run on a single CPU
#SBATCH --cpus-per-task=8                                      # Number of cores per task
#SBATCH --mem=40gb                                             # Job memory request
#SBATCH --time=07-00:00:00                                     # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err            # Standard error log

#SBATCH --mail-type=END,FAIL                                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                            # Where to send mail	

#Set output directory variable
OUTDIR="/scratch/ma95362/musse_MGA/merged"
REF="/scratch/ma95362/ET_1291/ref/gbk/genomic.gbk"

#Tell the program to make  the outdir folder
if [ ! -d $OUTDIR ] 
    then 
        mkdir -p $OUTDIR
fi

module load Bactopia/3.1.0
cd $OUTDIR
bactopia \
    --wf snippy \
    --reference $REF \
    --exclude $OUTDIR/bactopia-exclude.tsv \
    --bactopia $OUTDIR/MGA_paired_end_samples