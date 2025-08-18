#!/bin/bash
#SBATCH --job-name=Kraken2                                       # Job name
#SBATCH --partition=batch                                        # Partition (queue) name
#SBATCH --ntasks=1                                               # Run on a single CPU
#SBATCH --cpus-per-task=8                                        # Number of cores per task
#SBATCH --mem=40gb                                               # Job memory request
#SBATCH --time=05-00:00:00                                       # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err            # Standard error log

#SBATCH --mail-type=END,FAIL                                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                            # Where to send mail	

# Load Kraken2 module
module load Kraken2/2.1.3-gompi-2023a

# Paths
DB=/scratch/ma95362/k2_standard_16gb_20240605 # Path to Kraken2 database
INPUT=/scratch/ma95362/musse_MGA/merged/first_run_merged/all_reads/reads_for_Kraken2
OUTPUT=/scratch/ma95362/musse_MGA/merged/first_run_merged/all_reads/reads_for_Kraken2/kraken2_output

# Create output folder if it doesn't exist
mkdir -p $OUTPUT

cd $OUTPUT

# Loop over paired-end samples
for R1 in $INPUT/*_R1.fastq.gz
do
    SAMPLE=$(basename $R1 _R1.fastq.gz)
    R2=$INPUT/${SAMPLE}_R2.fastq.gz

    echo "Processing sample: $SAMPLE"

    kraken2 \
      --db $DB \
      --paired $R1 $R2 \
      --threads 8 \
      --report $OUTPUT/${SAMPLE}.report \
      --output $OUTPUT/${SAMPLE}.kraken
done

echo "All selected samples processed."
