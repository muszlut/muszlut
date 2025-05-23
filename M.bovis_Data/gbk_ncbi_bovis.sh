#!/bin/bash
#SBATCH --job-name=M.bovis_ref_gbk_ncbi       # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=8               # Number of cores per task
#SBATCH --mem=40gb                     # Job memory request
#SBATCH --time=03:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err             # Standard error log

#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu  # Where to send mail	

#Set output directory variable
OUTDIR="/scratch/ma95362/ETH_bovis_Sequence/bovis_REF"
#Set variable URL for genome website
URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/835/GCF_000195835.3_ASM19583v2/GCF_000195835.3_ASM19583v2_genomic.gbff.gz"

#Tell the program to make  the outdir folder
if [ ! -d $OUTDIR ] 
    then 
        mkdir -p $OUTDIR
fi
cd $OUTDIR
###Grab annotated m.canetti genome file from the NCBI URL variable then unzip it
curl -s $URL | gunzip -c > $OUTDIR/genomic.gbff
#module load NCBI-Datasets-CLI/16.4.4
#cd $OUTDIR
#ask merlin to run Staphy_ program 
#download datasets and unzip
#datasets download genome accession GCF_000195955.2 --include gbff
#unzip ncbi_dataset.zip
#addedddddd