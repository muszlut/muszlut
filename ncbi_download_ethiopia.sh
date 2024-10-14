#!/bin/bash
#SBATCH --job-name=mustutorial         # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=1gb                     # Job memory request
#SBATCH --time=02:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err             # Standard error log

#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu  # Where to send mail	

#Set output directory variable
OUTDIR="scratch/ma95362/ref"

#Tell the program to make  the outdir folder
if [ ! -d $OUTDIR ] 
    then 
        mkdir -p $OUTDIR
fi

 #Load modules
module load NCBI-Databasets-CLI/16.4.4

#move to working directory
cd $OUTDIR

#download datasets and unzip
datasets download genome accesssion GCF_000195955.2 --include cds, genome
unzip ncbi_download.zip

#ma95362@ss-sub3 ma95362$ pwd
#/scratch/ma95362
#ma95362@ss-sub3 ma95362$ cd /home/ma95362
#ma95362@ss-sub3 ~$ ls
#example  muszlut  newexample.txt  ondemand
#ma95362@ss-sub3 ~$ cd muszlut/
#ma95362@ss-sub3 muszlut$ git pull
#remote: Enumerating objects: 4, done.
#remote: Counting objects: 100% (4/4), done.
#remote: Compressing objects: 100% (2/2), done.
#remote: Total 3 (delta 1), reused 3 (delta 1), pack-reused 0 (from 0)
#Unpacking objects: 100% (3/3), 783 bytes | 195.00 KiB/s, done.
#From github.com:muszlut/muszlut
   dd7b9a5..170ab1f  main       -> origin/main
#Updating dd7b9a5..170ab1f
#Fast-forward
# ncbi_download_ethiopia.sh | 30 ++++++++++++++++++++++++++++++
# 1 file changed, 30 insertions(+)
# create mode 100644 ncbi_download_ethiopia.sh