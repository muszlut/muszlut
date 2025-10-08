#!/bin/bash
#SBATCH --job-name=Local_Pangenome
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out   # STDOUT log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err    # STDERR log
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
WORKDIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads"
BACTOPIA_DIR="${WORKDIR}/ETH_paired_end_samples"
FOFN="${WORKDIR}/all_samples.fofn"
OUTDIR="${WORKDIR}/pangenome_results"

mkdir -p "${OUTDIR}" "${WORKDIR}/logs"
cd "${OUTDIR}"

# -----------------------------
# Run Bactopia Pangenome Workflow
# -----------------------------
echo "[$(date)] Starting Bactopia Pangenome Workflow..."
bactopia \
    --wf pangenome \
    --bactopia "${BACTOPIA_DIR}" \
    --include "${FOFN}" \
    --outdir "${OUTDIR}" \
    --cpus ${SLURM_CPUS_PER_TASK} \
    --force \
    --skip_check_input \
    --species "Mycobacterium tuberculosis complex"

echo "[$(date)] ✅ Bactopia Pangenome completed successfully."

