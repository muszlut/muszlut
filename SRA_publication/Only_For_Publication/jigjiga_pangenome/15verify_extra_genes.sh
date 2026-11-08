#!/bin/bash
#SBATCH --job-name=mtb_verify_extra
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
# Put any additional candidate through the same verification as steps 12-14:
#   extract its sequence -> locate it in H37Rv -> read coverage from snippy's alignments
#   -> compare with the assembly call -> run the LMM on the read-based genotype.
#
# Written for group_862, a genome-wide significant EPTB hit (9.2e-5) that appears only in the
# panaroo-filter-pa matrix and has never been verified. Edit GENES to check others.
set -euo pipefail

GENES="${GENES:-group_862}"          # space-separated; override with GENES="a b" sbatch ...

module load Mamba/23.11.0-0
module load minimap2/2.28-GCCcore-13.2.0
source activate mtb_gwas

PROJ=/scratch/ma95362/Bio_project_publication
SCRIPTS=${SCRIPTS:-$PROJ/gwas/scripts}
PAN=$PROJ/PanPanaroo_merged_results/bactopia-runs/pangenome-20260927-192230/panaroo
PY=$PROJ/gwas/pyseer
GBKDIR=/scratch/ma95362/gbk/ncbi_dataset/data/GCF_000195955.2
OUT=$PROJ/gwas/verify_extra
MAN=$PROJ/gwas/depth_snippy/bam_manifest.tsv
mkdir -p "$OUT/depth"

[ -s "$MAN" ] || { echo "missing $MAN - run 13_depth_from_snippy.sh first"; exit 1; }

GENOME=$(ls -1 "$GBKDIR"/*.fna 2>/dev/null | head -1 || true)
[ -n "$GENOME" ] || GENOME=$PROJ/gwas/readlevel_wg/h37rv.fa
[ -s "$GENOME" ] || { echo "no H37Rv FASTA"; exit 1; }

echo "===== 1. extract $GENES from the pan-genome reference"
FA=$PAN/pan_genome_reference.fa
[ -s "$FA" ] || FA=$PAN/pan_genome_reference.fa.gz
[ -s "$FA" ] || { echo "no pan_genome_reference.fa(.gz) in $PAN"; exit 1; }
python3 - "$FA" "$OUT/extra.fa" $GENES <<'PY'
import gzip, sys
fa, out = sys.argv[1], sys.argv[2]
want = set(sys.argv[3:])
op = gzip.open if fa.endswith(".gz") else open
seqs, name, buf = {}, None, []
with op(fa, "rt") as f:
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            if name in want: seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else: buf.append(line)
    if name in want: seqs[name] = "".join(buf)
with open(out, "w") as o:
    for n in sorted(seqs):
        o.write(f">{n}\n")
        for i in range(0, len(seqs[n]), 70): o.write(seqs[n][i:i+70] + "\n")
        print(f"  {n}: {len(seqs[n])} bp")
for n in sorted(want - set(seqs)): print(f"  {n}: NOT FOUND in the pan-genome reference")
PY
[ -s "$OUT/extra.fa" ] || { echo "nothing extracted"; exit 1; }

echo; echo "===== 2. locate in H37Rv"
minimap2 -cx asm10 --secondary=yes -N 20 "$GENOME" "$OUT/extra.fa" 2>/dev/null \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$6,$8,$9,$10,$11,$12}' > "$OUT/hits.paf"
if [ ! -s "$OUT/hits.paf" ]; then
    echo "  NO alignment to H37Rv - this is an accessory gene absent from the reference."
    echo "  Coverage-at-coordinates cannot be used; it would need its own mapping. Stopping."
    exit 0
fi
printf "%-14s %-7s %-12s %-12s %-9s %-7s\n" gene len t_start t_end matches mapq
awk -F'\t' '{printf "%-14s %-7s %-12s %-12s %-9s %-7s\n", $1,$2,$4,$5,$6,$8}' "$OUT/hits.paf"
awk -F'\t' '{n[$1]++} END{for(g in n) if(n[g]>1) printf "  WARNING: %s hits %d H37Rv locations - not uniquely genotypable\n", g, n[g]}' "$OUT/hits.paf"

awk -F'\t' '{printf "%s\t%d\t%d\n", $1, $4+1, $5}' "$OUT/hits.paf" | sort -u > "$OUT/regions.tsv"
echo; echo "regions:"; cat "$OUT/regions.tsv"

echo; echo "===== 3. coverage from the existing snippy alignments"
: > "$OUT/depth_long.tsv"
printf "samples\tgene\tlength\tbreadth5_q20\tmeandepth_q20\tbreadth5_q0\tmeandepth_q0\n" > "$OUT/depth_wide.tsv"
n=0
while IFS=$'\t' read -r S BAM; do
    [ -s "$BAM" ] || continue
    CHR=$(samtools view -H "$BAM" | awk '$1=="@SQ"{sub("SN:","",$2); print $2; exit}')
    while IFS=$'\t' read -r G A B; do
        read -r b20 m20 <<<"$(samtools depth -aa -Q 20 -r "$CHR:$A-$B" "$BAM" 2>/dev/null \
            | awk '{n++; s+=$3; if($3>=5) c++} END{if(n) printf "%.4f %.2f",(c+0)/n,s/n; else print "0 0"}')"
        read -r b0 m0 <<<"$(samtools depth -aa -Q 0 -r "$CHR:$A-$B" "$BAM" 2>/dev/null \
            | awk '{n++; s+=$3; if($3>=5) c++} END{if(n) printf "%.4f %.2f",(c+0)/n,s/n; else print "0 0"}')"
        printf "%s\t%s\t%d\t%s\t%s\t%s\t%s\n" "$S" "$G" "$((B-A+1))" "$b20" "$m20" "$b0" "$m0" >> "$OUT/depth_wide.tsv"
        printf "%s\t%s\t%d\t%s\t%s\n" "$S" "$G" "$((B-A+1))" "$b20" "$m20" >> "$OUT/depth_long.tsv"
    done < "$OUT/regions.tsv"
    n=$((n+1)); [ $((n % 100)) -eq 0 ] && echo "  $n isolates done"
done < "$MAN"
echo "  $n isolates done"
awk -F'\t' -v d="$OUT/depth" '{print > (d"/"$1".tsv")}' "$OUT/depth_long.tsv"

echo; echo "===== 4. read-based vs assembly calls"
python3 "$SCRIPTS/readlevel_summary.py" "$OUT/depth" "$SCRIPTS/phenotypes.tsv" "$SCRIPTS/covariates.tsv" \
        "$PAN/gene_presence_absence.Rtab" "$PAN/struct_presence_absence.Rtab" \
        "$OUT/summary_extra.tsv" | tee "$OUT/report_extra.txt"

echo; echo "===== 5. MAPQ 20 vs MAPQ 0 (should agree for a single-copy locus)"
awk -F'\t' 'NR>1 {g=$2; n[g]++; if($4+0>=0.9) a[g]++; if($6+0>=0.9) b[g]++}
    END{for(g in n) printf "  %-14s Q20 %s/%s   Q0 %s/%s\n", g,(a[g]+0),n[g],(b[g]+0),n[g]}' \
    "$OUT/depth_wide.tsv" | sort

echo; echo "===== 6. LMM on the read-based genotype"
python3 "$SCRIPTS/make_readlevel_rtab.py" "$OUT/depth" "$SCRIPTS/covariates.tsv" "$OUT/extra"
if ! pyseer --phenotypes "$SCRIPTS/phenotypes.tsv" --phenotype-column EPTB --pres "$OUT/extra.Rtab" \
       --similarity "$PY/phylogeny_K.tsv" --lmm --min-af 0.01 --max-af 0.99 --cpu 2 \
       --output-patterns "$OUT/patterns_extra.txt" > "$OUT/pyseer_extra.tsv" 2> "$OUT/pyseer_extra.log"; then
    echo "  pyseer failed (a locus present in every isolate cannot be tested):"
    sed 's/^/      /' "$OUT/pyseer_extra.log"
else
    { column -t -s$'\t' "$OUT/pyseer_extra.tsv" 2>/dev/null || cat "$OUT/pyseer_extra.tsv"; } | cut -c1-150
fi
echo
echo "  reference points: p from the Panaroo Rtab and from the filtered matrix, and the"
echo "  genome-wide thresholds 2.58e-4 (Rtab, 194 patterns) and 1.29e-4 (filtered, 387)."
for v in $GENES; do
    pu=$(awk -F'\t' -v v="$v" '$1==v{print $4}' "$PY/pyseer_EPTB_genes.tsv" 2>/dev/null | head -1)
    pf=$(awk -F'\t' -v v="$v" '$1==v{print $4}' "$PROJ/gwas/filtered/pyseer_EPTB_genes_filtered.tsv" 2>/dev/null | head -1)
    printf "  %-14s p_Rtab=%-12s p_filtered=%-12s\n" "$v" "${pu:-NA}" "${pf:-NA}"
done

echo; echo "===== done ($(date)); results: $OUT"