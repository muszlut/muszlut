#!/bin/bash
#SBATCH --job-name=Sequence read download      # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=8                    # Run on a single CPU
#SBATCH --mem=64gb                     # Job memory request
#SBATCH --time=23:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err             # Standard error log

#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu  # Where to send mail	

#Set output directory variable
OUTDIR="/scratch/ma95362/M.bovis_Ethiopia/Sequence_reads"
#Set variable URL for genome website
URL=" https://outlookuga-my.sharepoint.com/:f:/r/personal/ma95362_uga_edu/Documents/64%20SH%20Survey%20Samples?csf=1&web=1&e=XvK9PU "

#Tell the program to make  the outdir folder
if [ ! -d $OUTDIR ] 
    then 
        mkdir -p $OUTDIR
fi
