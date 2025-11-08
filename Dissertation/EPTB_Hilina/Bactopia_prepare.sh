#!/bin/bash
#SBATCH --job-name=Bactopia_PRJNA1174701                        # Job name
#SBATCH --partition=batch                                       # Partition (queue) name
#SBATCH --ntasks=1                                              # Run on a single CPU
#SBATCH --cpus-per-task=4                                       # Number of cores per task
#SBATCH --mem=120gb                                              # Job memory request
#SBATCH --time=02-00:00:00                                      # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out            # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err             # Standard error log

#SBATCH --mail-type=END,FAIL                                    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                             # Where to send mail	

#Set output directory variable
OUTDIR="/scratch/ma95362/EPTB_Hilina/ETH_Bactopia_Prepare"

#Tell the program to make  the outdir folder
if [ ! -d $OUTDIR ] 
    then 
        mkdir -p $OUTDIR
fi

module load Bactopia/3.2.0

# Make output directory if it doesn't exist
mkdir -p $OUTDIR
cd $OUTDIR
bactopia prepare \
    --path /scratch/ma95362/EPTB_Hilina/reads \
    --species "Mycobacterium tuberculosis" \
    --genome-size 4410000 \
    > $OUTDIR/ETH_samples.txt
bactopia \
    --samples $OUTDIR/ETH_samples.txt \
    --coverage 100 \
    --outdir $OUTDIR/ETH_paired_end_samples \
    --max_cpus 4
bactopia summary \
    --bactopia-path $OUTDIR/ETH_paired_end_samples
