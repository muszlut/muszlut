#!/bin/bash
#SBATCH --job-name=Kraken_merge_py_script       # Job name
#SBATCH --partition=batch                       # Partition (queue) name
#SBATCH --ntasks=1                              # Run on a single CPU
#SBATCH --cpus-per-task=8                       # Number of cores per task
#SBATCH --mem=40gb                              # Job memory request
#SBATCH --time=03:00:00                         # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out   # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL                    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu            # Where to send mail

# Load modules
module load Python/3.11.5
module load Biopython/1.84-foss-2023b

# Set working and script directories
WORKDIR=/scratch/ma95362/ETH_bovis_Sequence/reads/64_SH_Sequence_data/raw/fixed_fastqs/fixed_fastqs
SCRIPT_DIR=/home/ma95362/muszlut/M.bovis_Data

# Move to working directory
cd $WORKDIR

# Run Python script
python $SCRIPT_DIR/Kraken_merge.py

