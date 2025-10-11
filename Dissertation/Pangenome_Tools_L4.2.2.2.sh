#!/bin/bash
#SBATCH --job-name=Pan_L4.2.2.2
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
OUTDIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/pangenome_tools_results/"

# Create output folder if not exist
mkdir -p $OUTDIR
cd $OUTDIR

# Run pangenome workflow using bactopia tools pangenome
#bactopia tools pangenome \
#    --samples /scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_samples.txt \
#    --exclude /scratch/ma95362/eth_national_analysis/all_fastq_reads/bactopia-exclude_final.tsv

#The error log from after runniging above command recommends to continue with the following command for pangenome analysis
bactopia -profile standard --bactopia bactopia --wf pangenome \
    --bactopia /scratch/ma95362/eth_national_analysis/all_fastq_reads/pangenome_tools_results/bactopia \
    --include /home/ma95362/muszlut/Dissertation/bactopia-exclude_except_L4.2.2.2.tsv