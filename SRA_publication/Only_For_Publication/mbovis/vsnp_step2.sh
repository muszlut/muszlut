#!/bin/bash
#SBATCH --job-name=vsnp3_step2
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/ma95362/vsnp3/logs/step2.%j.out
#SBATCH --error=/scratch/ma95362/vsnp3/logs/step2.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# =============================================================================
# vSNP3 step 2 - group analysis across the 38 samples from step 1.
#
# Assigns isolates to groups on defining SNPs -- this is what gives the
# M. bovis clonal complexes (Af1, Af2, Eu1, Eu2) that step 1 reported as
# "group file not provided" -- and produces annotated SNP tables, sorted
# matrices and phylogenetic trees.
#
# ---------------------------------------------------------------------------
# NOTES FROM EARLIER FAILURES
#   * step 2 consumes the ZERO-COVERAGE VCFs (*_zc.vcf), NOT the
#     *_filtered_hapall_annotated.vcf files. The zc VCF carries every
#     reference position, which is how vSNP distinguishes a true reference
#     call from missing data. Feeding it the annotated VCFs left
#     map_quality_dict empty and crashed with
#     "Length mismatch: Expected axis has 1 elements, new values have 2".
#   * vSNP3 truncates some sample names when writing the zc VCF
#     (097G_400 -> 097G_zc.vcf, 32_HK -> 32_zc.vcf). The VCFs are therefore
#     copied named after their step1 SAMPLE DIRECTORY so labels match
#     include.txt, TB-Profiler and the Snippy tree.
#   * step 2 ZIPS the input VCFs into vcf_starting_files.zip as it runs, so
#     a re-run finds no input and fails with
#     "'VCF_to_DF' object has no attribute 'chrom'". This script rebuilds
#     the VCF directory from scratch every time.
#   * 38 whole-genome (4.35 Mb) dataframes are held in memory at once and
#     vSNP3 spawns one worker per visible CPU. An 'interact' session (2 GB)
#     dies with BrokenProcessPool. Hence 180G / 8 cpus here.
#   * -m (metadata) is omitted: the USDA metadata sheet maps its own VCF
#     filenames to display names and will not match these samples.
#   * 'module load' must come BEFORE 'set -u' (Lmod's shell function
#     references unset variables internally).
#
# Run AFTER the step 1 array has finished:  sbatch vsnp3_step2.sh
# =============================================================================

module load vSNP3/3.31

set -euo pipefail

# ----------------------------- configuration ---------------------------------

BASE="/scratch/ma95362/vsnp3"
DEPS="/scratch/ma95362/vsnp3_test_dataset/vsnp_dependencies/Mycobacterium_AF2122"
STEP1="$BASE/step1"
STEP2="$BASE/step2"

DEFINING="$DEPS/Mbovis_define_filter.xlsx"
GBK="$DEPS/NC_002945v4.gbk"

# ------------------------------ preflight ------------------------------------

[ -s "$DEFINING" ] || { echo "ERROR: missing $DEFINING" >&2; exit 1; }
[ -s "$GBK" ]      || { echo "ERROR: missing $GBK" >&2; exit 1; }
[ -d "$STEP1" ]    || { echo "ERROR: missing $STEP1" >&2; exit 1; }

command -v vsnp3_step2.py >/dev/null 2>&1 \
    || { echo "ERROR: vsnp3_step2.py not on PATH after module load" >&2; exit 1; }

echo "Step 1 dir: $STEP1"
echo "Step 2 dir: $STEP2"
echo "Defining:   $DEFINING"
echo "Annotation: $GBK"
echo "CPUs:       ${SLURM_CPUS_PER_TASK:-?}   Memory: ${SLURM_MEM_PER_NODE:-?} MB"
echo "Started:    $(date)"
echo

# --------------------- rebuild the VCF input directory -----------------------
# Always from scratch: step 2 zips and removes its inputs as it runs, and a
# half-populated directory produces confusing downstream errors.

PREV="$STEP2.prev.$$"
[ -d "$STEP2" ] && mv "$STEP2" "$PREV" && echo "Previous step2 moved to $PREV"
mkdir -p "$STEP2"

missing=0
for d in "$STEP1"/*/; do
    s=$(basename "$d")
    v=$(find "$d" -name "*_zc.vcf" | head -1)
    if [ -n "$v" ]; then
        cp "$v" "$STEP2/${s}.vcf"
    else
        echo "WARNING: no *_zc.vcf found for $s" >&2
        missing=$((missing + 1))
    fi
done

n=$(ls -1 "$STEP2"/*.vcf 2>/dev/null | wc -l)
echo "Collected $n VCFs ($missing samples missing)"
[ "$n" -gt 1 ] || { echo "ERROR: need at least 2 VCFs, found $n" >&2; exit 1; }
echo

cd "$STEP2"

# --------------------------------- run step 2 --------------------------------
#   -a             table with ALL isolates as well as per-group tables
#   -s             defining SNPs -> group / clonal complex assignment
#   -b             GenBank annotation for the SNP tables
#   --show_groups  print group names in the SNP table
#   -html_tree     interactive HTML tree (also enables -dp, depth in tables)

vsnp3_step2.py \
    -a \
    -s "$DEFINING" \
    -b "$GBK" \
    --show_groups \
    -html_tree

echo
echo "Finished: $(date)"
echo
echo "Group directories (one per defining-SNP group / clonal complex):"
ls -d "$STEP2"/*/ 2>/dev/null || echo "  (none)"
echo
echo "Spreadsheets:"
ls -la "$STEP2"/*.xlsx 2>/dev/null || echo "  (none)"
echo
echo "Cross-check the groups against the TB-Profiler La1.2 / La1.3 / La2"
echo "calls and the SB numbers in $BASE/spoligotype_summary.csv"