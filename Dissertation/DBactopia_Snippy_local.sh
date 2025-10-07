#!/bin/bash
#SBATCH --job-name=Bactopia_Snippy_local_samples
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# -----------------------------
# Load environment modules
# -----------------------------
module load Bactopia/3.2.0-conda
module load Java/17.0.6


#Set output directory variable
OUTDIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/DBactopia_Snippy/Local_samples_only"
REF="/scratch/ma95362/gbk/ncbi_dataset/data/GCF_000195955.2/genomic.gbk"

# Create output folder if not exist
mkdir -p $OUTDIR
cd $OUTDIR
# Run snippy workflow using all samples
bactopia \
    --wf snippy \
    --reference $REF \
    --bactopia $OUTDIR/ETH_paired_end_samples \
    --include /scratch/ma95362/eth_national_analysis/all_fastq_reads/local_samples.fofn \
    --exclude $OUTDIR/ETH_paired_end_samples/bactopia-exclude.tsv