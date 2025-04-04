#!/bin/bash
#SBATCH --job-name=mtbvartools_first_attempt_Preprocessing     # Job name
#SBATCH --partition=batch                                      # Partition (queue) name
#SBATCH --ntasks=1                                             # Run on a single CPU
#SBATCH --cpus-per-task=8                                      # Number of cores per task
#SBATCH --mem=40gb                                             # Job memory request
#SBATCH --time=07-00:00:00                                     # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err            # Standard error log
#SBATCH --array=0-16                                           # This sets the job array for 17 isolates (0 to 16)

#SBATCH --mail-type=END,FAIL                                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                            # Where to send mail	

#Set output directory variable
OUTDIR="/scratch/ma95362/PRJNA823537_ET125/mtbvartool_first"
INPUT_DIR="/scratch/ma95362/PRJNA823537_ET125/ena-multiple-samples/fastqs"
###REF="/scratch/ma95362/ET_1291/ref/gbk/genomic.gbk"

#Tell the program to make  the outdir folder
if [ ! -d $OUTDIR ] 
    then 
        mkdir -p $OUTDIR
fi

#Load modules
module load Miniconda3/23.5.2-0

#activate conda envirnoment (should be version 1.0)
source activate mtbvartools

#move to working directory:
cd $OUTDIR

# List your 17 isolates or samples here
ISOLATES=("E22." "E23." "E38_M." "E61." "P105." "P122." "P124." "P14." "P147." "P209." "P42." "SRX1154627" "SRX14743237" "SRX2000723" "SRX7207533" "SRX8964326" "SRX9131782")

# Get the current isolate based on the SLURM_ARRAY_TASK_ID
ISOLATE="${ISOLATES[$SLURM_ARRAY_TASK_ID]}"

# Input files for each isolate (R1 and R2)
INPUT_FILE_R1="${INPUT_DIR}/${ISOLATE}_R1.fastq.gz"
INPUT_FILE_R2="${INPUT_DIR}/${ISOLATE}_R2.fastq.gz"

# Output directory for this isolate
OUTPUT_DIR_SAMPLE="${OUTPUT_DIR}/${ISOLATE}"

mkdir -p "$OUTPUT_DIR_SAMPLE"


# Run mtbvartools for the current isolate
mtbvartools --input $INPUT_FILE --output $OUTPUT_DIR_SAMPLE

echo "Analysis for ${ISOLATE} is complete."