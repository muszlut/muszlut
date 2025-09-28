#!/bin/bash
#SBATCH --job-name=bactopia_master
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8         # CPUs per Snippy job
#SBATCH --mem=32G                 # RAM per Snippy job
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/ma95362/logs/snippy_%A_%a.out
#SBATCH --error=/scratch/ma95362/logs/snippy_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu
#SBATCH --array=1-1398             # Adjust based on number of samples

# Load Bactopia
module load Bactopia/3.2.0

# -----------------------------
# Directories and files
# -----------------------------
READS_DIR=/scratch/ma95362/eth_national_analysis/all_fastq_reads
SNIPPY_DIR=/scratch/ma95362/eth_national_analysis/snippy
PANGENOME_DIR=/scratch/ma95362/eth_national_analysis/pangenome
SAMPLE_LIST=${READS_DIR}/sample_list.txt
REF="/scratch/ma95362/gbk/ncbi_dataset/data/GCF_000195955.2/genomic.gbk"

mkdir -p ${SNIPPY_DIR} ${PANGENOME_DIR}/roary ${PANGENOME_DIR}/panaroo

# -----------------------------
# Get sample for this array task
# -----------------------------
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_LIST})
echo "Processing sample: $SAMPLE"

# -----------------------------
# Step 1: Run Snippy workflow per sample
# -----------------------------
bactopia --reads "${READS_DIR}/${SAMPLE}_R1.fastq.gz" "${READS_DIR}/${SAMPLE}_R2.fastq.gz" \
         --wf snippy \
         --reference ${REF} \
         --genus Mycobacterium \
         --outdir ${SNIPPY_DIR}/${SAMPLE} \
         --cpus 8 \
         --mem 32 \
         --force

# -----------------------------
# Step 2: Pangenome workflows will be triggered automatically after all array jobs
# -----------------------------
if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
    # Submit dependent pangenome job only once (after all array jobs finish)
    jid=${SLURM_JOB_ID}
    sbatch --dependency=afterok:${jid} --export=ALL <<'EOF'
#!/bin/bash
#SBATCH --job-name=bactopia_pangenome
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/ma95362/logs/pangenome_%j.out
#SBATCH --error=/scratch/ma95362/logs/pangenome_%j.err

module load Bactopia/3.1.0

SNIPPY_DIR=/scratch/ma95362/eth_national_analysis/snippy
PANGENOME_DIR=/scratch/ma95362/eth_national_analysis/pangenome

# Run Roary
bactopia --wf pangenome
         --bactopia ${SNIPPY_DIR} \
         --outdir ${PANGENOME_DIR} \
         --cpus 32 \
         --mem 256 \
         --strict \
         --force

EOF
fi
echo "✅ Snippy job for ${SAMPLE} completed. Pangenome job will run after all samples are processed."