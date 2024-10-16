#!/bin/bash
#SBATCH --job-name=Bactopia_tutorial2        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=4               # Number of cores per task
#SBATCH --mem=40gb                     # Job memory request
#SBATCH --time=03:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err             # Standard error log

#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu  # Where to send mail	

#Set output directory variable
OUTDIR="/scratch/ma95362/fastqs"

#Tell the program to make  the outdir folder
if [ ! -d $OUTDIR ] 
    then 
        mkdir -p $OUTDIR
fi

module load Bactopia/3.1.0
cd $OUTDIR
bactopia prepare \
    --path $OUTDIR \
    --species "Staphylococcus aureus" \
    --genome-size 2800000
    > $OUTDIR/tutorial_samples.txt
bactopia \
    --samples $OUTDIR/tutorial_samples.txt \
    --coverage 100 \
    --outdir $OUTDIR/paired_end_samples\
    --max_cpus 4
bactopia summary \
    --bactopia-path $OUTDIR/paired_end_samples/