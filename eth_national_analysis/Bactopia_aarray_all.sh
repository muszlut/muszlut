#!/bin/bash
#SBATCH --job-name=bactopia_array_full
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8          # Enough for assembly and variant calling
#SBATCH --mem=32G                  # Reduced from 120G to allow more concurrent tasks
#SBATCH --time=05-00:00:00
#SBATCH --array=1-1399
#SBATCH --output=/scratch/ma95362/scratch/log.%A_%a.out
#SBATCH --error=/scratch/ma95362/scratch/log.%A_%a.err

module load Bactopia/3.2.0

# Path to samples.txt
SAMPLES=/scratch/ma95362/eth_national_analysis/bactopia_prepare/samples.txt

# Get the line for this task
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLES)

# Extract sample name and paired-end FASTQ files
SAMPLE=$(echo $LINE | cut -f1)
FQ1=$(echo $LINE | cut -f2)
FQ2=$(echo $LINE | cut -f3)

# Output directory for this sample
OUTDIR=/scratch/ma95362/eth_national_analysis/bactopia_results/$SAMPLE
mkdir -p $OUTDIR
cd $OUTDIR || { echo "Failed to cd into $OUTDIR"; exit 1; } 


# Skip already-processed samples
if [ -f "$OUTDIR/bactopia.done" ]; then
    echo "Sample $SAMPLE already processed. Skipping."
    exit 0
fi

# Run full Bactopia workflow on paired-end reads
bactopia --fq1 $FQ1 --fq2 $FQ2 \
         --outdir $OUTDIR \
         --species "Mycobacterium tuberculosis" \
         --cpus 8 \
         --genome-size 4400000

# Mark as done
touch $OUTDIR/bactopia.done
