#!/bin/bash
#SBATCH --job-name=vSNP_bovis
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/ma95362/vsnp_project/log.%j.out
#SBATCH --error=/scratch/ma95362/vsnp_project/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Load conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate vsnp3

# Base paths
BASE_DIR=/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples
STEP1_DIR=${BASE_DIR}/vsnp_project/step1_output
STEP2_DIR=${BASE_DIR}/vsnp_project/step2_output
REF_NAME=Mycobacterium_AF2122

mkdir -p $STEP1_DIR $STEP2_DIR

# Step 1: process each sample recursively
cd $STEP1_DIR
for R1 in $(find ${BASE_DIR} -type f -name "*_R1*.fastq.gz"); do
    R2=${R1/_R1/_R2}
    SAMPLE=$(basename ${R1} | cut -d"_" -f1)
    
    # Check that R2 exists
    if [[ ! -f "$R2" ]]; then
        echo "⚠️  Skipping $SAMPLE: R2 file not found"
        continue
    fi

    echo "▶ Running vSNP Step1 for sample: $SAMPLE"
    echo "   R1: $R1"
    echo "   R2: $R2"

    vsnp3_step1.py -r1 $R1 -r2 $R2 -t $REF_NAME -o ${STEP1_DIR}/${SAMPLE}_vsnp1
done

# Step 2: combine SNPs and build tree
cd $STEP2_DIR
vsnp3_step2.py -a -t $REF_NAME -i $STEP1_DIR

echo "✅ vSNP pipeline complete. Outputs in $STEP2_DIR"
