#!/bin/bash
#SBATCH --job-name=Bactopia_GACRC
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80gb
#SBATCH --time=07-00:00:00
#SBATCH --gres=lscratch:100
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# -------------------------------
# Variables
# -------------------------------
READS_DIR="/scratch/ma95362/clean_sequences_reads"
OUTDIR="/scratch/ma95362/SRAG_publication"
REF="/scratch/ma95362/gbk/ncbi_dataset/data/GCF_000195955.2/genomic.gbk"
# -------------------------------
# Load Bactopia
# -------------------------------
module purge
module load Bactopia/4.0.0-conda

# -------------------------------
# Use lscratch (CRITICAL FIX)
# -------------------------------
export BACTOPIA_CACHE=/lscratch/$SLURM_JOB_ID/bactopia-cache
mkdir -p "$BACTOPIA_CACHE"/{conda,datasets,singularity}

# Nextflow working directory (huge performance boost)
export NXF_WORK=/lscratch/$SLURM_JOB_ID/work
mkdir -p "$NXF_WORK"

# -------------------------------
# Create output directory
# -------------------------------
mkdir -p "$OUTDIR"

# -------------------------------
# Move to output directory
# -------------------------------
cd "$OUTDIR"

# -------------------------------
# Prepare samples
# -------------------------------
#bactopia prepare \
#    --path "$READS_DIR" \
#    --species "Mycobacterium tuberculosis" \
#    --genome-size 4410000 \
#    > "$OUTDIR/samples.txt"
#
# -------------------------------
# Run Bactopia
# -------------------------------
#bactopia \
#    --samples "$OUTDIR/samples.txt" \
#    --coverage 100 \
#    --outdir "$OUTDIR" \
#    --max_cpus 16 \
#    --max_memory 60.GB \
#    -resume
# Generate summary
# Clean bad version files
#find "$OUTDIR" -name "versions.yml" -delete
# Run summary
#bactopia summary --bactopia-path "$OUTDIR"

# Run Snippy workflow
# -------------------------------
bactopia \
    --wf snippy \
    --bactopia "$OUTDIR" \
    --reference "$REF" \
    --outdir "$OUTDIR/snippy1" \
    --max_cpus 16 \
    --max_memory 70.GB \
    -resume
 

