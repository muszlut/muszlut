#!/bin/bash
#SBATCH --job-name=Mtb_Eth763
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/ma95362/Eth_dataset_763/bactopia_run/logs/log.%j.out
#SBATCH --error=/scratch/ma95362/Eth_dataset_763/bactopia_run/logs/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# NOTE: create the logs folder BEFORE submitting, or Slurm cannot write the logs:
#   mkdir -p /scratch/ma95362/Eth_dataset_763/bactopia_run/logs

set -euo pipefail

# ============================================================
# 1. Environment
# ============================================================
module load Mamba/23.11.0-0
source activate bactopia4
export PYTHONNOUSERSITE=1

# Containers (recommended on HPC). Cache images once, reuse on every run.
module load Apptainer 2>/dev/null || module load Singularity 2>/dev/null || true
export NXF_SINGULARITY_CACHEDIR=/scratch/ma95362/singularity_cache
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

# ============================================================
# 2. Directories and settings
# ============================================================
FASTQ_DIR="/scratch/ma95362/Eth_dataset_763"          # where the FASTQs are
RUN_DIR="${FASTQ_DIR}/bactopia_run"
RESULT_DIR="${RUN_DIR}/results"
WORK_DIR="${RUN_DIR}/work"
SAMPLES="${RUN_DIR}/Eth763-samples.txt"

# Adjust these to match your file names, e.g.
#   ET001_R1.fastq.gz / ET001_R2.fastq.gz  -> FASTQ_EXT=.fastq.gz  (default patterns handle R1/R2)
#   ET001_1.fq.gz     / ET001_2.fq.gz      -> FASTQ_EXT=.fq.gz
FASTQ_EXT=".fastq.gz"
GENOME_SIZE=4411532        # H37Rv
SPECIES="Mycobacterium tuberculosis"

mkdir -p "$RESULT_DIR" "$WORK_DIR" "${RUN_DIR}/logs"
cd "$RUN_DIR"

echo "Bactopia version: $(bactopia --version)"

# ============================================================
# 3. Build the sample sheet (FOFN) from the local FASTQs
# ============================================================
echo "Preparing sample sheet from ${FASTQ_DIR} ..."

bactopia prepare \
    --path "$FASTQ_DIR" \
    --fastq-ext "$FASTQ_EXT" \
    --species "$SPECIES" \
    --genome-size "$GENOME_SIZE" \
    --recursive \
    > "$SAMPLES"

# ============================================================
# 4. Sanity checks on the sample sheet
# ============================================================
TOTAL=$(( $(wc -l < "$SAMPLES") - 1 ))     # minus header line
echo "=========================================="
echo "Samples in sheet: $TOTAL"
echo "Run types:"
tail -n +2 "$SAMPLES" | cut -f2 | sort | uniq -c
echo "=========================================="

if [[ "$TOTAL" -le 0 ]]; then
    echo "ERROR: no samples found. Check FASTQ_DIR and FASTQ_EXT."
    exit 1
fi

# Warn if anything was detected as single-end or needs merging (often a naming problem)
if tail -n +2 "$SAMPLES" | cut -f2 | grep -qvE '^paired-end$'; then
    echo "WARNING: some samples are not 'paired-end'. Inspect ${SAMPLES} before trusting results."
fi

# ============================================================
# 5. Run Bactopia
#    -resume lets you simply resubmit this script if the job
#    times out or fails; finished samples are not redone.
# ============================================================
echo "Starting Bactopia ..."

bactopia \
    --samples "$SAMPLES" \
    --coverage 100 \
    --max_cpus 32 \
    --max_memory 120.GB \
    --outdir "$RESULT_DIR" \
    -w "$WORK_DIR" \
    -profile singularity \
    -resume

echo "=========================================="
echo "Bactopia finished for $TOTAL samples"
echo "Results: $RESULT_DIR"
echo "=========================================="