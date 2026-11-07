#!/bin/bash
#SBATCH --job-name=vsnp3_step2
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/ma95362/vsnp3/logs/step2.%j.out
#SBATCH --error=/scratch/ma95362/vsnp3/logs/step2.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# =============================================================================
# vSNP3 step 2 - group analysis across the 38 samples from step 1.
#
# Groups isolates on defining SNPs -- this is what assigns the M. bovis
# clonal complexes (Af1, Af2, Eu1, Eu2) that step 1 could not provide
# ("Groups: group file not provided" in the per-sample stats) -- and
# produces annotated SNP tables, sorted matrices and phylogenetic trees.
#
# ---------------------------------------------------------------------------
# NOTES
#   * 'module load' must come BEFORE 'set -u': Lmod's module shell function
#     references unset variables internally and dies silently under -u.
#   * vsnp3_path_adder.py cannot register the reference (the module's
#     dependencies directory under /apps is read-only), so the three
#     dependency files are passed explicitly rather than via -t:
#         -s  defining SNPs  -> REQUIRED for clonal complex assignment
#         -m  metadata
#         -b  GenBank annotation
#   * step 1 writes VCFs into <sample>/alignment_NC_002945v4/, so the
#     collection step searches recursively. The zero-coverage VCFs
#     (*zc.vcf) are intermediates and must NOT be included.
#
# Run AFTER the step 1 array has finished:  sbatch vsnp3_step2.sh
# =============================================================================

module load vSNP3/3.31

set -euo pipefail

# ----------------------------- configuration ---------------------------------

BASE="/scratch/ma95362/vsnp3"
DEPS="/scratch/ma95362/vsnp3_test_dataset/vsnp_dependencies/Mycobacterium_AF2122"
STEP2="$BASE/step2"

DEFINING="$DEPS/Mbovis_define_filter.xlsx"
METADATA="$DEPS/Mbovis_metadata.xlsx"
GBK="$DEPS/NC_002945v4.gbk"

# ------------------------------ preflight ------------------------------------

for f in "$DEFINING" "$METADATA" "$GBK"; do
    [ -s "$f" ] || { echo "ERROR: missing dependency file: $f" >&2; exit 1; }
done

command -v vsnp3_step2.py >/dev/null 2>&1 \
    || { echo "ERROR: vsnp3_step2.py not on PATH after module load" >&2; exit 1; }

mkdir -p "$STEP2"
cd "$STEP2"

echo "Step 1 dir:  $BASE/step1"
echo "Step 2 dir:  $STEP2"
echo "Defining:    $DEFINING"
echo "Metadata:    $METADATA"
echo "Annotation:  $GBK"
echo "vSNP3:       $(command -v vsnp3_step2.py)"
echo "Started:     $(date)"
echo

# --------------------------- collect the VCFs --------------------------------

echo "Collecting VCFs from $BASE/step1 ..."
find "$BASE/step1" -name "*.vcf" ! -name "*zc.vcf" -exec cp -n {} . \;

n=$(ls -1 ./*.vcf 2>/dev/null | wc -l)
echo "Found $n VCF files"
[ "$n" -gt 1 ] || { echo "ERROR: need at least 2 VCFs, found $n" >&2; exit 1; }
echo

# --------------------------------- run step 2 --------------------------------
#   -a             build a table with ALL isolates as well as per-group tables
#   -s             defining SNPs -> clonal complex / group assignment
#   -m             metadata (sample renaming)
#   -b             GenBank annotation for the SNP tables
#   --show_groups  print group names in the SNP table
#   -html_tree     interactive HTML tree (also enables -dp, depth in tables)

vsnp3_step2.py \
    -a \
    -s "$DEFINING" \
    -m "$METADATA" \
    -b "$GBK" \
    --show_groups \
    -html_tree

echo
echo "Finished: $(date)"
echo
echo "Group directories (one per clonal complex / defining-SNP group):"
ls -d "$STEP2"/*/ 2>/dev/null || echo "  (none found)"
echo
echo "Spreadsheets produced:"
ls -la "$STEP2"/*.xlsx 2>/dev/null || echo "  (none found)"
echo
echo "The annotated SNP tables (*.xlsx) carry the group assignments and make"
echo "good supplementary material. Cross-check the groups against the"
echo "TB-Profiler La1.2 / La1.3 / La2 calls and the SB numbers in"
echo "$BASE/spoligotype_summary.csv"