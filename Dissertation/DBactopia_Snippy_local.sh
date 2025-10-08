#!/bin/bash
#SBATCH --job-name=Local_Snippy_
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# -----------------------------
# Load environment modules
# -----------------------------
module purge
module load Bactopia/3.2.0-conda
module load Java/17.0.6

# -----------------------------
# Define directories
# -----------------------------
OUTDIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads"
BACTOPIA_DIR="${OUTDIR}/ETH_paired_end_samples"
FOFN="${OUTDIR}/all_samples.fofn"
REF="/scratch/ma95362/gbk/ncbi_dataset/data/GCF_000195955.2/genomic.gbk"
RESULTS="${OUTDIR}/snippy_results"

mkdir -p "${RESULTS}"
cd "${RESULTS}"

# -----------------------------
# Run Bactopia Snippy Workflow (resumable)
# -----------------------------
echo "[$(date)] Starting Bactopia Snippy Workflow (resumable mode)..."

bactopia \
    --wf snippy \
    --reference "${REF}" \
    --bactopia "${BACTOPIA_DIR}" \
    --include "${FOFN}" \
    --outdir "${RESULTS}" \
    --cpus ${SLURM_CPUS_PER_TASK} \
    --skip_check_input \
    --species "Mycobacterium tuberculosis complex" \
    -resume

echo "[$(date)] ✅ Bactopia Snippy completed (or resumed) successfully."
