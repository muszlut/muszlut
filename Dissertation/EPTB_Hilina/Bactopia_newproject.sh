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
# Step 0: Set Output Directory
# -------------------------------
OUTDIR="/scratch/ma95362/EPTB_Hilina/new_project"
mkdir -p $OUTDIR
cd $OUTDIR

echo "---- Step 1: Loading Modules ----"
module purge
module load Bactopia/3.2.0-conda
module load EDirect/20.5.20231006-GCCcore-12.3.0    # needed for esearch/efetch

# -------------------------------
# Step 2: Fetch SRR accessions
# -------------------------------
echo "---- Step 2: Fetching SRR accessions from BioProject PRJNA1174701 ----"

esearch -db sra -query PRJNA1174701 \
    | efetch -format runinfo \
    | cut -d',' -f1 \
    | grep SRR \
    > srr_list.txt

echo "Found $(wc -l < srr_list.txt) SRR accessions."
echo "Saved to: $OUTDIR/srr_list.txt"

# safety check
if [[ ! -s srr_list.txt ]]; then
    echo "ERROR: No SRR IDs extracted. Exiting."
    exit 1
fi

# -------------------------------
# Step 3: Run Bactopia using SRR IDs
# -------------------------------
echo "---- Step 3: Running Bactopia on SRR IDs ----"

bactopia \
    --accessions srr_list.txt \
    --outdir $OUTDIR \
    --force \
    --max_cpus 16 \
    --max_memory "120 GB" \
    --max_time "72h"

echo "---- Bactopia RUN COMPLETE ----"

