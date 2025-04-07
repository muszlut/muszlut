#!/bin/bash
#SBATCH --job-name=Mtbvartools_job                           # Job name
#SBATCH --partition=batch                                    # Partition (queue) name
#SBATCH --ntasks=1                                           # Run on a single CPU
#SBATCH --cpus-per-task=8                                    # Number of cores per task
#SBATCH --mem=40gb                                           # Job memory request
#SBATCH --time=03:00:00                                      # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out         # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err          # Standard error log

#SBATCH --mail-type=END,FAIL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                          # Where to send mail

module load Python/3.10.4  # Load the Python module
source ~/mtbvartools_env/bin/activate  # Activate your virtual environment (if applicable)

python -c "import mtbvartools; print('mtbvartools is working!')"