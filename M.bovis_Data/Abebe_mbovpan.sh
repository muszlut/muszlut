#!/bin/bash
#SBATCH --job-name=Abebe_mbovpan_pipeline                    # Job name
#SBATCH --partition=batch                                    # Partition (queue) name
#SBATCH --ntasks=1                                           # Run on a single CPU
#SBATCH --cpus-per-task=32                                   # Number of cores per task
#SBATCH --mem=80gb                                           # Job memory request (Increased from 40GB)
#SBATCH --time=05-00:00:00                                   # Time limit d-hh:mm:ss
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out         # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err          # Standard error log
#SBATCH --mail-type=END,FAIL                                 # Mail notifications
#SBATCH --mail-user=ma95362@uga.edu                          # Your email address

# IMPORTANT: Ensure there are NO leading spaces or invisible characters
# before any #SBATCH directives on these lines.

OUTDIR="/scratch/ma95362/ETH_bovis_Sequence/Abebe_samples/mbovpan_result"

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Activate conda environment
source ~/.bashrc
conda activate mbovpan-env

# Navigate to mbovpan working directory
cd /scratch/ma95362/mbovpan

# Run the pipeline in full mode (spoligotyping + SNP + pangenome + virulence gene profiling)
# Nextflow will now use the process-specific resource settings from nextflow.config
nextflow run main.nf \
  --input /scratch/ma95362/ETH_bovis_Sequence/mbovpanfastq_reads \
  --output "$OUTDIR" \
  --run all \
  --threads 32 \
  --qual 20 \
  --depth 25 \
  --mapq 40