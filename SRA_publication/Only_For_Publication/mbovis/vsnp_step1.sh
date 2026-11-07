#!/bin/bash
#SBATCH --job-name=vsnp3_step1
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=1-38%10
#SBATCH --output=/scratch/ma95362/vsnp3/logs/step1.%A_%a.out
#SBATCH --error=/scratch/ma95362/vsnp3/logs/step1.%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# =============================================================================
# vSNP3 step 1 - per-sample alignment, variant calling and spoligotyping
# against M. bovis AF2122/97 (NC_002945.4).
#
# Runs as a SLURM array, one task per sample from include.txt (38 samples,
# 10 concurrent). vSNP3 has no scheduler integration, so parallelism is at
# the job level.
#
# ---------------------------------------------------------------------------
# NOTES FROM EARLIER FAILURES
#   * 'module load' must come BEFORE 'set -u': Lmod's module shell function
#     references unset variables internally and dies silently under -u.
#   * The reference files are now named explicitly. The previous version used
#     'ls "$DEPS"/*.fasta "$DEPS"/*.fna | head -1'; with no .fna present, ls
#     exited 2, 'set -o pipefail' propagated it, and 'set -e' killed the
#     script before any output -- exit code 2, empty logs, 2 seconds.
#   * The module's dependencies directory under /apps is read-only, so
#     vsnp3_path_adder.py cannot register a reference type. The FASTA and
#     GBK are therefore passed explicitly with -f and -b instead of -t.
#
# Run:  sbatch vsnp3_step1.sh
# Test: sbatch --array=1-1 vsnp3_step1.sh
# =============================================================================

module load vSNP3/3.31

set -euo pipefail

# ----------------------------- configuration ---------------------------------

READS="/scratch/ma95362/Bio_project_publicatio/mbovis_seq"
BASE="/scratch/ma95362/vsnp3"
INCLUDE="/scratch/ma95362/Bio_project_publication/Mbovis/include.txt"
DEPS="/scratch/ma95362/vsnp3_test_dataset/vsnp_dependencies/Mycobacterium_AF2122"

REF_FASTA="$DEPS/NC_002945v4.fasta"
REF_GBK="$DEPS/NC_002945v4.gbk"

mkdir -p "$BASE/logs" "$BASE/step1"

# ------------------------------ preflight ------------------------------------

[ -s "$INCLUDE"   ] || { echo "ERROR: $INCLUDE missing or empty" >&2; exit 1; }
[ -s "$REF_FASTA" ] || { echo "ERROR: $REF_FASTA missing" >&2; exit 1; }
[ -s "$REF_GBK"   ] || { echo "ERROR: $REF_GBK missing" >&2; exit 1; }

command -v vsnp3_step1.py >/dev/null 2>&1 \
    || { echo "ERROR: vsnp3_step1.py not on PATH after module load" >&2; exit 1; }

# --------------------------- pick this task's sample -------------------------

: "${SLURM_ARRAY_TASK_ID:?must be run as a SLURM array job}"

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INCLUDE" | tr -d '[:space:]')
[ -n "$SAMPLE" ] || { echo "ERROR: no sample at line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }

R1="$READS/${SAMPLE}_R1.fastq.gz"
R2="$READS/${SAMPLE}_R2.fastq.gz"
[ -s "$R1" ] || { echo "ERROR: $R1 missing" >&2; exit 1; }
[ -s "$R2" ] || { echo "ERROR: $R2 missing" >&2; exit 1; }

WORK="$BASE/step1/$SAMPLE"
mkdir -p "$WORK"
cd "$WORK"

echo "Sample:     $SAMPLE  (array task ${SLURM_ARRAY_TASK_ID})"
echo "Reference:  $REF_FASTA"
echo "Annotation: $REF_GBK"
echo "Workdir:    $WORK"
echo "vSNP3:      $(command -v vsnp3_step1.py)"
echo "Started:    $(date)"
echo

# vSNP3 expects the FASTQs in the working directory
ln -sf "$R1" .
ln -sf "$R2" .

# ------------------------------- run step 1 ----------------------------------
# -spoligo is required: spoligotyping for TB complex isolates is NOT run by
# default in vSNP3 3.x.

vsnp3_step1.py \
    -n "$SAMPLE" \
    -r1 "${SAMPLE}_R1.fastq.gz" \
    -r2 "${SAMPLE}_R2.fastq.gz" \
    -f "$REF_FASTA" \
    -b "$REF_GBK" \
    -spoligo

echo
echo "Finished: $(date)"
echo
ls -la "$WORK"