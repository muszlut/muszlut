#!/bin/bash
#SBATCH --job-name=bactopia_array_all       # Job name
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=64              # CPUs per job
#SBATCH --mem=230G               # highmem per job
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err             # Standard error log

#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu  # Where to send mail	
#SBATCH --array=1-1398           # one job per sample

module load Bactopia/3.1.0

# Directories
READS_DIR=/scratch/ma95362/eth_national_analysis/all_fastq_reads
PROJECT_DIR=/scratch/ma95362/eth_national_analysis/project
PANGENOME_DIR=/scratch/ma95362/eth_national_analysis/pangenome
SNIPPY_DIR=/scratch/ma95362/eth_national_analysis/snippy
SAMPLE_LIST=/scratch/ma95362/eth_national_analysis/all_fastq_reads/sample_list.txt
REF="/scratch/ma95362/gbk/ncbi_dataset/data/GCF_000195955.2/genomic.gbk"

# Create output directories if they don't exist
mkdir -p $PROJECT_DIR
mkdir -p $PANGENOME_DIR
mkdir -p $SNIPPY_DIR

# Get sample name
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_LIST})

# Step 1: Run Bactopia per sample (QC + assembly + annotation)
bactopia --reads "${READS_DIR}/${SAMPLE}_R1.fastq.gz" "${READS_DIR}/${SAMPLE}_R2.fastq.gz" \
         --genus Mycobacterium \
         --outdir ${PROJECT_DIR}/${SAMPLE} \
         --cpus 64

# Step 2: Run pangenome only once after all array jobs are done
# Use SLURM job dependencies to trigger pangenome after array completion
# Submit pangenome job separately:
# sbatch --dependency=afterok:<ARRAY_JOB_ID> pangenome_highmem.slurm
