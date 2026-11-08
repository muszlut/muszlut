#!/bin/bash
#SBATCH --job-name=mtb_masked_tree
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu
#
# Masked core-SNP phylogeny for 257 M. tuberculosis isolates
#  1. mask PE/PPE, PE_PGRS, IS/transposase and phage CDS (from H37Rv GenBank)
#  2. apply the mask to the snippy-core whole-genome alignment, drop Reference
#  3. keep variable sites with complete ACGT calls (snp-sites -c)
#  4. IQ-TREE: ModelFinder (+ASC), 1000 ultrafast bootstraps + 1000 SH-aLRT
set -euo pipefail

module load Mamba/23.11.0-0
source activate mtb_gwas

PROJ=/scratch/ma95362/Bio_project_publication
SNIPPY=$PROJ/snippy_results/bactopia-runs/snippy-20260927-161232
GBK=/scratch/ma95362/gbk/ncbi_dataset/data/GCF_000195955.2/genomic.gbk   # H37Rv NC_000962.3
SCRIPTS=$PROJ/gwas/scripts
OUT=$PROJ/gwas/tree
mkdir -p "$OUT"; cd "$OUT"

[ -s "$GBK" ] || { echo "H37Rv GenBank not found at $GBK"; exit 1; }
[ -s "$SNIPPY/core-snp-clean.full.aln.gz" ] || { echo "snippy alignment not found in $SNIPPY"; exit 1; }

echo "== 1. mask"
python3 "$SCRIPTS/make_mask_from_gbk.py" "$GBK" > H37Rv_mask.bed
wc -l H37Rv_mask.bed

echo "== 2. apply mask (whole-genome alignment, 4,411,532 bp)"
python3 "$SCRIPTS/mask_alignment.py" "$SNIPPY/core-snp-clean.full.aln.gz" H37Rv_mask.bed \
        core-snp.masked_PEPPE.full.aln --drop-reference

echo "== 3. variable sites"
snp-sites -c -o core-snp.masked_PEPPE.snps.aln core-snp.masked_PEPPE.full.aln
seqkit stats core-snp.masked_PEPPE.snps.aln

echo "== 4. IQ-TREE"
IQTREE=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree)
echo "using $IQTREE ($($IQTREE --version 2>&1 | head -1))"
$IQTREE -s core-snp.masked_PEPPE.snps.aln -m MFP+ASC -B 1000 -alrt 1000 -T 16 \
        --prefix mtb257_masked --seed 2026

gzip -f core-snp.masked_PEPPE.full.aln
echo "== done: $OUT/mtb257_masked.treefile"
grep -E "Best-fit model|Number of sequences|Number of constant|parsimony-informative|Log-likelihood of the tree" mtb257_masked.iqtree