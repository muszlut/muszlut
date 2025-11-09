#!/bin/bash
#SBATCH --job-name=Bactopia_Mtb
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120gb
#SBATCH --time=03-00:00:00
#SBATCH --output=/scratch/ma95362/log.%j.out
#SBATCH --error=/scratch/ma95362/log.%j.err

#SBATCH --mail-type=END,FAIL                                   
#SBATCH --mail-user=ma95362@uga.edu                         

module load Bactopia/3.1.0

OUTDIR="/scratch/ma95362/EPTB_Hilina/ETH_Bactopia_Prepare"
READS="/scratch/ma95362/EPTB_Hilina/reads"

# Step 1: Prepare
bactopia prepare \
  --path $READS \
  --species "Mycobacterium tuberculosis" \
  --genome-size 4410000 \
  > $OUTDIR/ETH_samples.txt

# Step 2: Run analysis
bactopia \
  --samples $OUTDIR/ETH_samples.txt \
  --coverage 100 \
  --outdir $OUTDIR/ETH_paired_end_samples \
  --max_cpus 8

# Step 3: Generate summary
bactopia summary \
  --bactopia-path $OUTDIR/ETH_paired_end_samples

