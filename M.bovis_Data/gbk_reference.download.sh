#!/bin/bash
#SBATCH --job-name=gbk_ncbi       # Job name
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
OUTDIR="/home/ma95362/.conda/envs/mtbvartools/share/snpeff-5.2-1/data/AF2122_97"
#Set variable URL for genome website
URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/835/GCF_000195835.3_ASM19583v2/GCF_000195835.3_ASM19583v2_genomic.gff.gz"

#Tell the program to make  the outdir folder
if [ ! -d $OUTDIR ] 
    then 
        mkdir -p $OUTDIR
fi

###Grab annotated mtb H37Rv genome file from the NCBI URL variable then unzip it
curl -s $URL | gunzip -c > $OUTDIR/genomic.gff