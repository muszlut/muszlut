#!/bin/bash
#SBATCH --job-name=Abebe_mbovpan_pipeline
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=80gb
#SBATCH --time=05-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

OUTDIR="/scratch/ma95362/ETH_bovis_Sequence/Abebe_samples/mbovpan_result"

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Activate conda environment
source ~/.bashrc
conda activate mbovpan-env

# Navigate to mbovpan working directory
cd /scratch/ma95362/mbovpan

# Run the pipeline in full mode (spoligotyping + SNP + pangenome + virulence gene profilingg)
nextflow run main.nf \
  --input /scratch/ma95362/ETH_bovis_Sequence/mbovpanfastq_reads \
  --output "$OUTDIR" \
  --run all \
  --threads 32 \
  --qual 20 \
  --depth 25 \
  --mapq 40
