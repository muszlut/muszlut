#!/bin/bash
#SBATCH --job-name=mtb_followup
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/ma95362/Bio_project_publication/logs/log.%j.out
#SBATCH --error=/scratch/ma95362/Bio_project_publication/logs/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu
#
# Follow-up on the GWAS hits: T3-ETH marker SNP / gene carriers, cas2-cas1-csm6
# arrangement, EPTB candidate carriers, candidate sequences for BLAST, and
# resistance-gene SNPs in the tested SNP set. Results -> gwas/followup/
set -euo pipefail

module load Mamba/23.11.0-0
source activate mtb_gwas

PROJ=/scratch/ma95362/Bio_project_publication
SCRIPTS=${SCRIPTS:-$PROJ/gwas/scripts}
PAN=$PROJ/PanPanaroo_merged_results/bactopia-runs/pangenome-20260927-192230/panaroo
GBK=/scratch/ma95362/gbk/ncbi_dataset/data/GCF_000195955.2/genomic.gbk
SNPS=$PROJ/gwas/T3ETH/snps_matrix.tsv
OUT=$PROJ/gwas/followup
mkdir -p "$OUT"; cd "$OUT"

for f in "$SCRIPTS/followup.py" "$SCRIPTS/sample_key.tsv" "$SNPS" "$PAN/gene_presence_absence.Rtab" \
         "$PAN/struct_presence_absence.Rtab" "$GBK"; do
    [ -s "$f" ] || { echo "missing input: $f"; exit 1; }
done

python3 "$SCRIPTS/followup.py" "$SCRIPTS/sample_key.tsv" "$SNPS" "$PAN" "$GBK" "$OUT" | tee followup_report.txt

echo; echo "===== done ($(date)); report: $OUT/followup_report.txt; sequences: $OUT/followup_sequences.fasta"