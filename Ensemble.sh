#!/bin/bash
#SBATCH --job-name=Ensemble download ethiopia        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=1gb                     # Job memory request
#SBATCH --time=02:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=cdlog.%j.err             # Standard error log

#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu  # Where to send mail	

#Set output directory variable
OUTDIR="/scratch/ma95362/ref"
#Set variable URL for genome website
URL=" https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-59/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.59.gff3.gz "

#Tell the program to make  the outdir folder
if [ ! -d $OUTDIR ] 
    then 
        mkdir -p $OUTDIR
fi

###Grab annotated genome file from the URL variable the unzip it and put it in a file called ecoli_MG1655.gff
curl -s $URL | gunzip -c > $OUTDIR/ecoli_MG1655.gff
