#!/bin/bash
#SBATCH --job-name=Musse_Spoligotyping                     # Job name
#SBATCH --partition=batch                                  # Partition (queue) name
#SBATCH --ntasks=1                                         # Run on a single CPU
#SBATCH --cpus-per-task=8                                  # Number of cores per task
#SBATCH --mem=40gb                                         # Job memory request
#SBATCH --time=03-00:00:00                                 # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out       # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err        # Standard error log
#SBATCH --mail-type=END,FAIL                               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                        # Where to send mail

# Variables
OUTDIR="/scratch/ma95362/musse_MGA/fastqs"
FASTQ_DIR="/scratch/ma95362/musse_MGA/fastqs/MGA_paired_end_samples"  # Directory containing FASTQ files
FOFN='/scratch/ma95362/musse_MGA/fastqs/MGA_samples.txt'                 # File containing sample names
EXCLUDE_FILE='/scratch/ma95362/musse_MGA/fastqs/bactopia-exclude.tsv'    # File containing samples to exclude

# Create output directory if it doesn't exist
if [ ! -d $OUTDIR ]; then
    mkdir -p $OUTDIR
fi

# Load the Bactopia module
module load Bactopia/3.1.0

# Create a set of samples to exclude
EXCLUDES=$(cat $EXCLUDE_FILE)

# Run the Bactopia spoligotyping workflow for each sample
while read SAMPLE; do 
    if echo "$EXCLUDES" | grep -qw "$SAMPLE"; then
        echo "Excluding sample $SAMPLE"
    else
        bactopia \
            -profile singularity \
            --wf tbprofiler \
            --input /work/fdqlab/Ethiopia_wgs_mtb_2024/first_run/${SAMPLE}.bam \
            --output ${OUTDIR}/${SAMPLE} \
            --spoligotype
    fi
done < $FOFN
