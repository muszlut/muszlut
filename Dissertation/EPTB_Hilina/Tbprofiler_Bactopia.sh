#!/bin/bash
#SBATCH --job-name=gbk_ncbi       # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=8               # Number of cores per task
#SBATCH --mem=40gb                     # Job memory request
#SBATCH --time=06:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err             # Standard error log

#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu  # Where to send mail	
#Set output directory variableklkjl
OUTDIR="/scratch/ma95362/EPTB_Hilina/Bactopia_Run"


#Tell the program to make  the outdir folder
if [ ! -d $OUTDIR ] 
    then 
        mkdir -p $OUTDIR
fi
module load Bactopia/3.2.0-conda
module load EDirect/20.5.20231006-GCCcore-12.3.0
cd $OUTDIR

bactopia \
    -profile singularity \
    --wf tbprofiler \
    --bactopia $OUTDIR/ETH_full_analysis