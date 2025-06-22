#!/bin/bash
#SBATCH --job-name=Abebe_Tbprofiler_Lineages_Spoligo         # Job name
#SBATCH --partition=batch                                    # Partition (queue) name
#SBATCH --ntasks=1                                           # Run on a single CPU
#SBATCH --cpus-per-task=8                                    # Number of cores per task
#SBATCH --mem=40gb                                           # Job memory request
#SBATCH --time=05-00:00:00                                   # Time limit d-hh:mm:ss
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out         # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err          # Standard error log
#SBATCH --array=0-39                                         # Update this based on the number of BAMs - 1
#SBATCH --mail-type=END,FAIL                                 # Mail notifications
#SBATCH --mail-user=ma95362@uga.edu                          # Your email address

# Load conda environment
source ~/.bashrc
conda activate tb-profiler-env

# Set directories
BAM_DIR=/scratch/ma95362/ETH_bovis_Sequence/Abebe_all_samples/all_bams
OUT_DIR=/scratch/ma95362/ETH_bovis_Sequence/Abebe_all_samples/tbprofiler_output

# Get list of BAM files
BAM_FILES=($BAM_DIR/*.bam)
BAM_FILE=${BAM_FILES[$SLURM_ARRAY_TASK_ID]}
SAMPLE=$(basename "$BAM_FILE" .bam)

#move to working directory
cd $OUTDIR

# Run TBProfiler
tb-profiler profile --bam "$BAM_FILE" --prefix "$OUT_DIR/$SAMPLE"
