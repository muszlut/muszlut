#!/bin/bash
#SBATCH --job-name=Bactopia_Pangenome        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=8               # Number of cores per task
#SBATCH --mem=40gb                     # Job memory request
#SBATCH --time=05-00:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err             # Standard error log

#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu  # Where to send mail	

# Set directories
FASTQ_DIR="/scratch/ma95362/PRJNA1056148_bactopia/ena-multiple-samples/"
OUT_DIR="/scratch/ma95362/Pangenome_ETH/"

#Tell the program to make  the outdir folder
if [ ! -d $OUT_DIR ] 
    then 
        mkdir -p $OUT_DIR
fi

# Load Bactopia and necessary tools
module load Bactopia

# Run Bactopia for pangenome analysis
bactopia --fastqs $FASTQ_DIR --datasets datasets/ --workdir $OUT_DIR \
  --tools pangenome \
  --cpus 16


echo "Pangenome analysis started"
