#!/bin/bash
#SBATCH --job-name=Musse_WGS_BacSpoligo                        # Job name
#SBATCH --partition=batch                                      # Partition (queue) name
#SBATCH --ntasks=1                                             # Run on a single CPU
#SBATCH --cpus-per-task=8                                      # Number of cores per task
#SBATCH --mem=40gb                                             # Job memory request
#SBATCH --time=03-00:00:00                                     # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err            # Standard error log

#SBATCH --mail-type=END,FAIL                                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                            # Where to send mail 

# Define directories
#WORKDIR="/scratch/ma95362/musse_MGA/fastqs/MGA_paired_end_samples"   # Directory where Bactopia is run
OUTDIR="/scratch/ma95362/musse_MGA/fastqs/MGA_paired_end_samples"    # Directory for Bactopia results
READSDIR="/work/fdqlab/Ethiopia_wgs_mtb_2024/first_run"      # Directory containing FASTQ files
DBDIR="/scratch/ma95362/musse_MGA/fastqs/MGA_paired_end_samples"   # Path to the Bactopia database

# Set output directory variable
OUTDIR="/scratch/ma95362/musse_MGA/fastqs/MGA_paired_end_samples"

# Tell the program to make the outdir folder
if [ ! -d $OUTDIR ]; then 
    mkdir -p $OUTDIR
fi

# Load necessary modules (if needed on Sapelo2)
module load bactopia

# Activate conda environment for Bactopia
#source activate bactopia


# Run Bactopia pipeline with spoligotyping enabled
bactopia \
    --R1 "$READSDIR/sample_1.fastq.gz" \
    --R2 "$READSDIR/sample_2.fastq.gz" \
    --datasets "$DBDIR" \
    --outdir "$OUTDIR" \
    --tools tb-profiler \
    --conda-env "$(conda info --base)/envs/bactopia-env"

echo "Spoligotyping analysis with Bactopia completed."
