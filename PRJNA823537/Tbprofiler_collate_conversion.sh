#!/bin/bash
#SBATCH --job-name=Tbprofiler_jason_collate_conversion       # Job name
#SBATCH --partition=batch                                    # Partition (queue) name
#SBATCH --ntasks=1                                            # Run on a single CPU
#SBATCH --cpus-per-task=8                                     # Number of cores per task
#SBATCH --mem=40gb                                             # Job memory request
#SBATCH --time=03:00:00                                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out              # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err             # Standard error log

#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu  # Where to send mail	
#Set output directory variable
OUTDIR="/scratch/ma95362/PRJNA823537_ET125/Tbprofiler_module/results"
SCRIPT="/home/ma95362/muszlut/PRJNA823537"

#Load modules 
module load Biopython/1.84-foss-2023b


#move to working diectory
cd $OUTDIR

# Run your script
python ~/Tbprofiler_collate_from_new_module.py --dir results/