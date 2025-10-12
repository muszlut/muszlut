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

# -----------------------------
# Base directories
BASE_DIR=/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples
STEP1_DIR=${BASE_DIR}/vsnp_project/step1_output
STEP2_DIR=${BASE_DIR}/vsnp_project/step2_output
mkdir -p $STEP1_DIR $STEP2_DIR

# Reference
REF_DIR=/home/ma95362/vsnp3_test_dataset/vsnp_dependencies/Mycobacterium_AF2122/Fasta
REF_FASTA=${REF_DIR}/NC_002945v4.fasta
if [[ ! -f "$REF_FASTA" ]]; then
    echo "❌ Reference FASTA not found: $REF_FASTA"
    exit 1
fi
REF_NAME=Mycobacterium_AF2122   # matches folder in vsnp_dependencies

# -----------------------------
# Step 1: Process each sample recursively
cd $STEP1_DIR
for R1 in $(find ${BASE_DIR} -type f -name "*_R1*.fastq.gz"); do
    R2=${R1/_R1/_R2}
    SAMPLE=$(basename ${R1} | cut -d"_" -f1)

    # Skip if R2 is missing
    if [[ ! -f "$R2" ]]; then
        echo "⚠️  Skipping $SAMPLE: R2 not found"
        continue
    fi

    echo "▶ Running vSNP Step1 for sample: $SAMPLE"
    echo "   R1: $R1"
    echo "   R2: $R2"
    echo "   Reference: $REF_FASTA"

    vsnp3_step1.py -r1 $R1 -r2 $R2 -f $REF_FASTA -o ${STEP1_DIR}/${SAMPLE}_vsnp1
done

# -----------------------------
# Step 2: Combine SNPs and build tree
cd $STEP2_DIR
# Use either reference folder (-t) or FASTA (-f)
vsnp3_step2.py -a -f $REF_FASTA

echo "✅ vSNP pipeline complete. Outputs in $STEP2_DIR"
