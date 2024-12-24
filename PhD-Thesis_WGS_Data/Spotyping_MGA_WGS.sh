#!/bin/bash
#SBATCH --job-name=Spotyping_Multiple                 # Job name
#SBATCH --partition=batch                             # Partition (queue) name
#SBATCH --ntasks=1                                    # Run on a single CPU
#SBATCH --cpus-per-task=8                             # Number of cores per task
#SBATCH --mem=40gb                                    # Job memory request
#SBATCH --time=03-00:00:00                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out  # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err   # Standard error log
#SBATCH --mail-type=END,FAIL                          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                   # Where to send mail

# Set output directory variable
OUTDIR="/scratch/ma95362/musse_MGA/spotyping_results"
FOFN='/scratch/ma95362/musse_MGA/fastqs/MGA_samples.txt'
EXCLUDE_FILE='/scratch/ma95362/musse_MGA/fastqs/bactopia-exclude.tsv'

# Tell the program to make the outdir folder if it doesn't exist
if [ ! -d $OUTDIR ]; then 
    mkdir -p $OUTDIR
fi

# Load Miniconda or the specific environment containing SpoTyping
module load Miniconda3/23.5.2-0

# Activate the environment with SpoTyping installed
source activate spotyping-env

# Create a set of samples to exclude
EXCLUDES=$(cat $EXCLUDE_FILE)

# Run SpoTyping for each sample, excluding specified ones
while read SAMPLE; do
    if echo "$EXCLUDES" | grep -qw "$SAMPLE"; then
        echo "Excluding sample $SAMPLE"
    else
        spotype.py -i /work/fdqlab/Ethiopia_wgs_mtb_2024/first_run/${SAMPLE}_R1.fastq.gz -r /work/fdqlab/Ethiopia_wgs_mtb_2024/first_run/${SAMPLE}_R2.fastq.gz -o ${OUTDIR}/${SAMPLE}_spoligotyping_results.txt
    fi
done < $FOFN
