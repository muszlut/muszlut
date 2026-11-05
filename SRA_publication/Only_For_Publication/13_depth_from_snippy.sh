#!/bin/bash
#SBATCH --job-name=mtb_depth_snippy
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/ma95362/Bio_project_publication/logs/log.%j.out
#SBATCH --error=/scratch/ma95362/Bio_project_publication/logs/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu
#
# THE CORRECT MEASUREMENT, replacing steps 08 and 11.
#
# Step 12 established that all five candidates occur exactly ONCE in H37Rv, full length, at
# MAPQ 60. So reads can be assigned unambiguously and the right question is simply: how deeply
# is each candidate's locus covered in each isolate's whole-genome alignment?
#
# Those alignments already exist - snippy mapped every isolate to H37Rv - so nothing is
# remapped here. Both MAPQ 0 and MAPQ 20 are reported; with single-copy loci the two should
# agree closely, and a large gap would itself indicate a mapping problem.
set -euo pipefail

module load Mamba/23.11.0-0
source activate mtb_gwas

PROJ=/scratch/ma95362/Bio_project_publication
SCRIPTS=${SCRIPTS:-$PROJ/gwas/scripts}
HITS=$PROJ/gwas/locate/candidate_hits.tsv
OUT=$PROJ/gwas/depth_snippy
mkdir -p "$OUT"

[ -s "$HITS" ] || { echo "missing $HITS - run 12_locate_candidates.sh first"; exit 1; }

# ---- 1. one BAM per isolate ----
# BAMs live in nextflow work directories and the same isolate may appear more than once
# (snippy plus other steps). Prefer a path mentioning snippy/variants; otherwise take the
# first. Only isolates listed in covariates.tsv are kept.
MAN=$OUT/bam_manifest.tsv
if [ ! -s "$MAN" ]; then
    echo "building BAM manifest ..."
    find "$PROJ" -name '*.bam' 2>/dev/null | grep -v -E 'readlevel|readlevel_wg' > "$OUT/all_bams.txt" || true
    : > "$MAN"
    while read -r s; do
        best=$(grep -F "/$s.bam" "$OUT/all_bams.txt" | grep -i -m1 -E 'snippy|variant' || true)
        [ -n "$best" ] || best=$(grep -F -m1 "/$s.bam" "$OUT/all_bams.txt" || true)
        [ -n "$best" ] && printf "%s\t%s\n" "$s" "$best" >> "$MAN"
    done < <(tail -n +2 "$SCRIPTS/covariates.tsv" | cut -f1)
    echo "isolates with a BAM: $(wc -l < "$MAN") of 257"
fi
head -2 "$MAN"
miss=$(( 257 - $(wc -l < "$MAN") ))
[ "$miss" -eq 0 ] || echo "WARNING: $miss isolates have no BAM and will be absent from the table"

# ---- 2. regions, converted from the PAF (0-based, end-exclusive) to 1-based inclusive ----
REG=$OUT/regions.tsv
awk -F'\t' 'NR>1 {printf "%s\t%d\t%d\n", $1, $6+1, $7}' "$HITS" | sort -u > "$REG"
echo; echo "regions:"; cat "$REG"

# ---- 3. depth per isolate per locus, at MAPQ 0 and 20 ----
echo; echo "reading coverage from the existing alignments ..."
: > "$OUT/depth_long.tsv"
printf "samples\tgene\tlength\tbreadth5_q20\tmeandepth_q20\tbreadth5_q0\tmeandepth_q0\n" > "$OUT/depth_wide.tsv"
n=0
while IFS=$'\t' read -r S BAM; do
    [ -s "$BAM" ] || continue
    [ -s "$BAM.bai" ] || samtools index "$BAM" 2>/dev/null || true
    CHR=$(samtools view -H "$BAM" | awk '$1=="@SQ"{sub("SN:","",$2); print $2; exit}')
    while IFS=$'\t' read -r G A B; do
        read -r br20 md20 <<<"$(samtools depth -aa -Q 20 -r "$CHR:$A-$B" "$BAM" 2>/dev/null \
            | awk '{n++; s+=$3; if($3>=5) c++} END{if(n) printf "%.4f %.2f", (c+0)/n, s/n; else print "0 0"}')"
        read -r br0 md0 <<<"$(samtools depth -aa -Q 0 -r "$CHR:$A-$B" "$BAM" 2>/dev/null \
            | awk '{n++; s+=$3; if($3>=5) c++} END{if(n) printf "%.4f %.2f", (c+0)/n, s/n; else print "0 0"}')"
        printf "%s\t%s\t%d\t%s\t%s\t%s\t%s\n" "$S" "$G" "$((B-A+1))" "$br20" "$md20" "$br0" "$md0" \
            >> "$OUT/depth_wide.tsv"
        # long form in the layout readlevel_summary.py expects (Q20)
        printf "%s\t%s\t%d\t%s\t%s\n" "$S" "$G" "$((B-A+1))" "$br20" "$md20" >> "$OUT/depth_long.tsv"
    done < "$REG"
    n=$((n+1)); [ $((n % 50)) -eq 0 ] && echo "  $n isolates done"
done < "$MAN"
echo "  $n isolates done"

# per-isolate files, so readlevel_summary.py can be reused unchanged
mkdir -p "$OUT/depth"
awk -F'\t' -v d="$OUT/depth" '{print > (d"/"$1".tsv")}' "$OUT/depth_long.tsv"

echo; echo "===== summary (presence = >=90% of the locus covered at >=5x, MAPQ>=20)"
PAN=$PROJ/PanPanaroo_merged_results/bactopia-runs/pangenome-20260927-192230/panaroo
python3 "$SCRIPTS/readlevel_summary.py" "$OUT/depth" "$SCRIPTS/phenotypes.tsv" "$SCRIPTS/covariates.tsv" \
        "$PAN/gene_presence_absence.Rtab" "$PAN/struct_presence_absence.Rtab" \
        "$OUT/summary_snippy.tsv" | tee "$OUT/report_snippy.txt"

echo; echo "===== MAPQ 20 vs MAPQ 0 (single-copy loci: these should agree)"
awk -F'\t' 'NR>1 {g=$2; n[g]++; if($4+0>=0.9) a[g]++; if($6+0>=0.9) b[g]++}
    END{printf "  %-14s %-16s %-16s\n","gene","present Q20","present Q0";
        for(g in n) printf "  %-14s %-16s %-16s\n", g, (a[g]+0)"/"n[g], (b[g]+0)"/"n[g]}' \
    "$OUT/depth_wide.tsv" | sort

echo; echo "===== done ($(date)); tables: $OUT"