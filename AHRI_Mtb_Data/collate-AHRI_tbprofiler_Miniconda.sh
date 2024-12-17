#!/bin/bash
#SBATCH --job-name=collate_Tbprofiler       # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=2               # Number of cores per task
#SBATCH --mem=10gb                     # Job memory request
#SBATCH --time=01:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err             # Standard error log

#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu  # Where to send mail	
#Set output directory variable
OUTDIR="/scratch/ma95362/PRJNA1155881_AHRI_bactopia/ena-multiple-samples"
FOFN='/scratch/ma95362/PRJNA1155881_AHRI_bactopia/ena-multiple-samples/bactopia-runs/snippy-20241108-111044'


#Tell the program to make  the outdir folder
if [ ! -d $OUTDIR ] 
    then 
        mkdir -p $OUTDIR
fi

#Load modules
module load Miniconda3/23.5.2-0

#activate conda envirnoment (should be version 6.3.0)

source activate tb-profiler-env

#move to working directory:
cd $OUTDIR

#samples just needs to be a list of sample names. No path is required.
tb-profiler collate --samples $FOFN/genomic.samples.txt --dir $OUTDIR/*/tools/tbprofiler --itol

