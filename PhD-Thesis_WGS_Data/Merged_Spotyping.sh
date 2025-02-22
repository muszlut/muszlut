#!/bin/bash
#SBATCH --job-name=Musse_Merged_SpoTyping                      # Job name
#SBATCH --partition=batch                                      # Partition (queue) name
#SBATCH --ntasks=1                                             # Run on a single CPU
#SBATCH --cpus-per-task=8                                      # Number of cores per task
#SBATCH --mem=40gb                                             # Job memory request
#SBATCH --time=05-00:00:00                                     # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err            # Standard error log

#SBATCH --mail-type=END,FAIL                                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                            # Where to send mail	

#Set output directory variable
OUTDIR="/scratch/ma95362/musse_MGA/merged/SpoTyping_results"

#Tell the program to make  the outdir folder
if [ ! -d $OUTDIR ] 
    then 
        mkdir -p $OUTDIR
fi

# Load required modules on Sapelo2
module load Bactopia/3.1.0
module load Singularity/3.11.4-GCCcore-11.3.0
module load conda

# Step 1: Run Bactopia to generate assemblies
bactopia \
  --samples samplesheet.csv \
  --path /work/fdqlab/merged_reads/Ethiopia_wgs_mtb_2024/first_run \
  --outdir ./bactopia_results \
  --wf myco \                           # Use "myco" workflow for Mycobacterium
  --genome_size 4.4M \                  # M. tuberculosis genome size
  --species "Mycobacterium tuberculosis" \
  -profile singularity \                # Use Singularity for containers
 
# Step 2: Set up SpoTyping in a Conda environment
conda create -n spotyping -y python=3.8 blast=2.13.0
source activate spotyping

# Clone SpoTyping repository
git clone https://github.com/xiaeryu/SpoTyping.git


for ASSEMBLY in bactopia_results/assemblies/*.fna.gz; do
  # Decompress assembly (required by SpoTyping)
  BASENAME=$(basename ${ASSEMBLY} .fna.gz)
  gunzip -c ${ASSEMBLY} > ${BASENAME}.fna

  # Run SpoTyping
  python SpoTyping/SpoTyping.py \
    -s ${BASENAME}.fna \
    -o spotyping_results/${BASENAME}

  # Clean up decompressed files
  rm ${BASENAME}.fna
done

# Step 4: Aggregate results (optional)
echo "Sample,Spoligotype,SIT" > spotyping_results/spotyping_summary.csv
for RESULT in spotyping_results/*/spotyping.txt; do
  SAMPLE=$(basename $(dirname ${RESULT}))
  SPOLIGOTYPE=$(awk '/Binary/ {print $3}' ${RESULT})
  SIT=$(awk '/SIT/ {print $3}' ${RESULT})
  echo "${SAMPLE},${SPOLIGOTYPE},${SIT}" >> spotyping_results/spotyping_summary.csv
done