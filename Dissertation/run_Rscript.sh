#!/bin/bash
#SBATCH --job-name=Rscript_T3_ETH        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=4              # Number of cores per task
#SBATCH --mem=40gb                     # Job memory request
#SBATCH --time=03:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err             # Standard error log

#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu  # Where to send mail	

#Set output directory variable
OUTDIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/pangenome_tools_results/bactopia/bactopia-runs/pangenome_of_L4.2.2.2.2/panaroo/filtered_output/pyseer_T3ETH_out"
#Location of R script
TIDY="/home/ma95362/muszlut/Dissertation/T3_ETH.R"


#Move to output directory
cd $OUTDIR

#Activate Renvironmemt previously created
source activate r-tidyverse

#Run R script:here
R --no-save < $TIDY
