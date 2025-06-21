#!/bin/bash
#SBATCH --job-name=Abebe_Mbovis_Snippy                       # Job name
#SBATCH --partition=batch                                    # Partition (queue) name
#SBATCH --ntasks=1                                           # Run on a single CPU
#SBATCH --cpus-per-task=8                                    # Number of cores per task
#SBATCH --mem=40gb                                           # Job memory request
#SBATCH --time=05-00:00:00                                   # Time limit d-hh:mm:ss
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out         # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err          # Standard error log
#SBATCH --mail-type=END,FAIL                                 # Mail notifications
#SBATCH --mail-user=ma95362@uga.edu                          # Your email address

# Set variables
OUTDIR="/scratch/ma95362/ETH_bovis_Sequence/Abebe_all_samples"
REF="/scratch/ma95362/ETH_bovis_Sequence/bovis_REF/genomic.gbk"
SNIPPY_OUT="$OUTDIR/snippy_results"

# Create output directory if it doesn’t exist
mkdir -p "$SNIPPY_OUT"

# Load Bactopia
module load Bactopia/3.1.0

#move to working directory
cd $OUTDIR

# Run snippy workflow
bactopia \
  --wf snippy \
  --reference "$REF" \
  --bactopia "$OUTDIR" \
  --outdir "$SNIPPY_OUT" \
  --max_cpus 8
