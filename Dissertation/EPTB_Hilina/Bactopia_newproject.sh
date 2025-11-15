#!/bin/bash
#SBATCH --job-name=Bactopia_Mtb
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120gb
#SBATCH --time=03-00:00:00
#SBATCH --output=/scratch/ma95362/EPTB_Hilina/new_project_logs/log.%j.out
#SBATCH --error=/scratch/ma95362/EPTB_Hilina/new_project_logs/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# -------------------------------
# SETUP OUTPUT DIRECTORY
# -------------------------------
OUTDIR="/scratch/ma95362/EPTB_Hilina/new_project"
mkdir -p $OUTDIR
cd $OUTDIR

echo "---- Step 1: Loading Modules ----"
module load Bactopia/3.2.0-conda
module load EDirect      # Needed to fetch SRR accessions

# -------------------------------
# STEP 1: FETCH SRR ACCESSIONS FROM BIOPROJECT
# -------------------------------
echo "---- Step 2: Fetching SRR accessions from BioProject PRJNA1174701 ----"

esearch -db sra -query PRJNA1174701 \
    | efetch -format runinfo \
    | cut -d ',' -f 1 \
    | grep SRR \
    > srr_list.txt

NUM_SRR=$(wc -l < srr_list.txt)

echo "Found $NUM_SRR SRR accessions."
echo "Saved to: ${OUTDIR}/srr_list.txt"

if [ "$NUM_SRR" -eq 0 ]; then
    echo "ERROR: No SRR entries found! Check the BioProject or network connection."
    exit 1
fi

# -------------------------------
# STEP 2: RUN BACTOPIA USING SRR IDs
# -------------------------------
echo "---- Step 3: Running Bactopia on SRR IDs ----"

bactopia \
    --accessions srr_list.txt \
    --outdir $OUTDIR \
    --cpus 16 \
    --genome-size 4400000 \
    --min-coverage 20 \
    --force

echo "---- Bactopia finished successfully ----"
