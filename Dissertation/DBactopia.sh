#!/bin/bash
#SBATCH --job-name=DBactopia_Run
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=23
#SBATCH --mem=240G
#SBATCH --time=07-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out    # STDOUT log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err     # STDERR log

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

#Set output directory variable
OUTDIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads"

#Tell the program to make  the outdir folder
if [ ! -d $OUTDIR ] 
    then 
        mkdir -p $OUTDIR
fi


# Load Bactopia
module load Bactopia/3.2.0-conda

cd $OUTDIR

bactopia prepare \
    --path $OUTDIR \
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
bactopia plot \
    --bactopia-path $OUTDIR/ETH_paired_end_samples \
    --outdir $OUTDIR/ETH_paired_end_samples/plots