#!/bin/bash
#SBATCH --job-name=Ref_setup                                   # Job name
#SBATCH --partition=batch                                      # Partition (queue) name
#SBATCH --ntasks=1                                             # Run on a single CPU
#SBATCH --cpus-per-task=8                                      # Number of cores per task
#SBATCH --mem=8gb                                             # Job memory request
#SBATCH --time=05-00:00:00                                     # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err            # Standard error log

#SBATCH --mail-type=END,FAIL                                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                            # Where to send mail	

OUTDIR="/scratch/ma95362/new_magma_project/reference"

#Tell the program to make  the outdir folder
if [ ! -d $OUTDIR ] 
    then 
        mkdir -p $OUTDIR
fi

cd $OUTDIR


# Download reference
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.gbff.gz

gunzip GCF_000195955.2_ASM19595v2_genomic.fna.gz
gunzip GCF_000195955.2_ASM19595v2_genomic.gbff.gz

mv GCF_000195955.2_ASM19595v2_genomic.fna reference.fasta
mv GCF_000195955.2_ASM19595v2_genomic.gbff reference.gbk

module load BWA
module load GATK
module load SAMTOOLS

bwa index reference.fasta
gatk CreateSequenceDictionary -R reference.fasta
samtools faidx reference.fasta