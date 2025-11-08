#!/bin/bash
#SBATCH --job-name=Bactopia_PRJNA1174701                        # Job name
#SBATCH --partition=batch                                       # Partition (queue) name
#SBATCH --ntasks=1                                              # Run on a single CPU
#SBATCH --cpus-per-task=8                                       # Number of cores per task
#SBATCH --mem=40gb                                              # Job memory request
#SBATCH --time=02-00:00:00                                      # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out            # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err             # Standard error log

#SBATCH --mail-type=END,FAIL                                    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                             # Where to send mail	

#Set output directory variable
OUTDIR="/scratch/ma95362/EPTB_Hilina/"

#Tell the program to make  the outdir folder
if [ ! -d $OUTDIR ] 
    then 
        mkdir -p $OUTDIR
fi

module load Bactopia/3.2.0-conda
module load Java/17.0.6

# Make output directory if it doesn't exist
mkdir -p $OUTDIR
cd $OUTDIR


bactopia prepare \
    --path /scratch/ma95362/eth_national_analysis/all_fastq_reads \
    --species "Mycobacterium tuberculosis" \
    --genome-size 4400000 \
    > /scratch/ma95362/eth_national_analysis/bactopia_prepare/samples.txt
