#!/bin/bash
#SBATCH --job-name=Bactopia_PRJNA1155881        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=8               # Number of cores per task
#SBATCH --mem=40gb                     # Job memory request
#SBATCH --time=02-00:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err             # Standard error log

#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu  # Where to send mail	

#Set output directory variable
OUTDIR="/scratch/ma95362/PRJNA1155881_AHRI_bactopia"

#Tell the program to make  the outdir folder
if [ ! -d $OUTDIR ] 
    then 
        mkdir -p $OUTDIRcd
fi

module load Bactopia/3.1.0
cd $OUTDIR
bactopia summary --bactopia-path $OUTDIR/ena-multiple-samples/
#bactopia search \
#    --query PRJNA1155881
#bactopia \
#    --accessions $OUTDIR/bactopia-accessions.txt \
#    --coverage 100 \
#    --outdir $OUTDIR/ena-multiple-samples \
#    --max_cpus 8