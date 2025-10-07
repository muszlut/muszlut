#!/bin/bash
#SBATCH --job-name=Bactopia_Pangenome_local
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


# Set directories
OUTDIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/"
#FOFN="/scratch/ma95362/eth_national_analysis/all_fastq_reads/local_samples.fofn"
# Create output folder if not exist
mkdir -p $OUTDIR
cd $OUTDIR

# Run pangenome workflow using FOFN
bactopia \
    --wf pangenome \
    --bactopia $OUTDIR/ETH_paired_end_samples \
    --include /scratch/ma95362/eth_national_analysis/all_fastq_reads/local_samples.fofn 
