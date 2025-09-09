#!/bin/bash
#SBATCH --job-name=Bactopia_query_p                     # Job name
#SBATCH --partition=batch                                # Partition (queue)
#SBATCH --ntasks=1                                       # Single task
#SBATCH --cpus-per-task=8                                # CPUs per task
#SBATCH --mem=40gb                                       # Memory
#SBATCH --time=05-00:00:00                               # Time limit (HH:MM:SS)
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out    # STDOUT log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err     # STDERR log
#SBATCH --mail-type=END,FAIL                             # Mail events
#SBATCH --mail-user=ma95362@uga.edu                      # Your email

#----------------------------
# Set output directory
#----------------------------
OUTDIR="/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads"
mkdir -p $OUTDIR

#----------------------------
# Load Bactopia module
#----------------------------
module load Bactopia/3.2.0
cd $OUTDIR
#----------------------------
# Run Bactopia query

bactopia search \
    --query PRJEB32192
bactopia \
    --accessions $OUTDIR/bactopia-accessions.txt \
    --coverage 100 \
    --outdir $OUTDIR/M.bovis_paired_end_samples \
    --max_cpus 8 \
    --skip_qc false \
    --skip_kraken2 false    # Explicitly run Kraken2 for species ID
#----------------------------
# Generate summary report
#----------------------------
bactopia summary \
    --bactopia-path $OUTDIR/M.bovis_paired_end_samples/