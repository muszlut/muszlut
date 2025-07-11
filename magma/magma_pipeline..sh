#!/bin/bash
#SBATCH --job-name=Mbovis_MAGMA               # job name
#SBATCH --partition=batch                    # partition (queue) name
#SBATCH --ntasks=1                           # run on a single task
#SBATCH --cpus-per-task=8                    # number of CPU cores per task
#SBATCH --mem=48G                            # memory requested
#SBATCH --time=05-00:00:00                   # time limit (days-hh:mm:ss)
#SBATCH --output=logs/magma_%j.out           # standard output log
#SBATCH --error=logs/magma_%j.err            # standard error log
#SBATCH --mail-type=END,FAIL                 # email notifications
#SBATCH --mail-user=ma95362@uga.edu          # your email address

# Define working/project directory
WORKDIR="/scratch/ma95362/ETH_bovis_Sequence/magma_output/work"
NXF_HOME="/scratch/ma95362/ETH_bovis_Sequence/magma_output/.nextflow"

# Create necessary directories
mkdir -p logs
mkdir -p "$WORKDIR"
mkdir -p "$NXF_HOME"

# Load Nextflow
module load Nextflow

# Export environment variables so Nextflow uses these directories
export NXF_WORK="$WORKDIR"
export NXF_HOME="$NXF_HOME"

# Run MAGMA
nextflow run 'https://github.com/TORCH-Consortium/MAGMA' \
    -profile conda_local,server \
    -r v1.1.1 \
    -params-file my_parameters_1.yml
