#!/bin/bash
#SBATCH --job-name=mtb_filtered
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/ma95362/Bio_project_publication/logs/log.%j.out
#SBATCH --error=/scratch/ma95362/Bio_project_publication/logs/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu
#
# Rerun the gene-based association tests on the Panaroo-FILTERED presence/absence matrix.
#
# Panaroo flags clusters that are likely annotation/assembly artefacts:
#   _stop  potential pseudogenes
#   _len   genes of unusual length for their cluster
#   ;      fragmented genes (multiple sequence IDs)
# panaroo-filter-pa removes them. This is Panaroo's own answer to the problem we
# demonstrated by read mapping in step 08, so running it (a) confirms that finding with the
# tool's native filter and (b) pre-empts the reviewer question "why not panaroo-filter-pa?".
#
# Also inspects the gene-graph neighbourhood of the cas locus and the EPTB candidates.
set -euo pipefail

module load Mamba/23.11.0-0
source activate mtb_gwas

PROJ=/scratch/ma95362/Bio_project_publication
SCRIPTS=${SCRIPTS:-$PROJ/gwas/scripts}
PAN=$PROJ/PanPanaroo_merged_results/bactopia-runs/pangenome-20260927-192230/panaroo
PY=$PROJ/gwas/pyseer
OUT=$PROJ/gwas/filtered
mkdir -p "$OUT"; cd "$OUT"

for f in "$PAN/gene_presence_absence.csv" "$SCRIPTS/phenotypes.tsv" "$PY/phylogeny_K.tsv"; do
    [ -s "$f" ] || { echo "missing input: $f"; exit 1; }
done

echo "===== 1. how many clusters are flagged in the unfiltered matrix"
tot=$(( $(wc -l < "$PAN/gene_presence_absence.csv") - 1 ))
awk -F',' 'NR>1{n++; if($1 ~ /_stop/) s++; if($1 ~ /_len/) l++}
           END{printf "  clusters: %d   _stop: %d   _len: %d\n", n, s+0, l+0}' \
    "$PAN/gene_presence_absence.csv"
echo "  rows containing a fragmented-gene separator (;): $(grep -c ';' "$PAN/gene_presence_absence.csv" || true)"

echo; echo "===== 2. are our candidates flagged?"
for g in group_1905 group_1965 group_589 group_1946 group_2137; do
    hit=$(grep -m1 -E "^\"?${g}(_stop|_len)?\"?," "$PAN/gene_presence_absence.csv" || true)
    if [ -n "$hit" ]; then
        name=$(echo "$hit" | cut -d',' -f1)
        frag=$(echo "$hit" | grep -c ';' || true)
        printf "  %-12s cluster name: %-22s fragmented-entries: %s\n" "$g" "$name" "$frag"
    else
        printf "  %-12s not found by exact name (may be merged into another cluster)\n" "$g"
    fi
done

echo; echo "===== 3. panaroo-filter-pa (remove pseudogenes, length outliers, fragmented)"
# panaroo-filter-pa builds its output path as <out_dir> + <the input path it was given>, so
# an absolute -i produces a nonsensical destination. Copy the matrix in and use a bare name.
cp -f "$PAN/gene_presence_absence.csv" "$OUT/gene_presence_absence.csv"
( cd "$OUT" && panaroo-filter-pa -i gene_presence_absence.csv -o ./ --type pseudo,length,frag )
ls -1 "$OUT"/gene_presence_absence*.Rtab 2>/dev/null || true
RTAB=$(ls -1 "$OUT"/*_filt_*.Rtab 2>/dev/null | head -1)
[ -n "$RTAB" ] || RTAB=$(ls -1 "$OUT"/*filt*.Rtab "$OUT"/*.Rtab 2>/dev/null | head -1)
[ -n "$RTAB" ] || { echo "could not find the filtered Rtab; check the panaroo-filter-pa output above"; exit 1; }
echo "filtered matrix: $RTAB"
echo "  clusters before: $tot   after: $(( $(wc -l < "$RTAB") - 1 ))"

echo; echo "===== 4. rerun the gene tests on the filtered matrix"
npat () { sort -u "$1" | wc -l; }
run () {   # $1 phenotype column  $2 phenotype file  $3 min-af
    local tag="$1_genes_filtered"
    if ! pyseer --phenotypes "$2" --phenotype-column "$1" --pres "$RTAB" \
           --similarity "$PY/phylogeny_K.tsv" --lmm \
           --min-af "$3" --max-af 0.95 --cpu 8 \
           --output-patterns "patterns_${tag}.txt" > "pyseer_${tag}.tsv" 2> "pyseer_${tag}.log"; then
        echo "--- $1: pyseer FAILED:"; sed 's/^/      /' "pyseer_${tag}.log"; return 0
    fi
    local n thr nsig
    n=$(npat "patterns_${tag}.txt"); thr=$(awk -v n="$n" 'BEGIN{printf "%.3g", 0.05/n}')
    nsig=$(awk -F'\t' -v t="$thr" 'NR>1 && $4!="" && $4<t' "pyseer_${tag}.tsv" | wc -l)
    printf -- "--- %-22s patterns=%-6s threshold=%-10s significant=%s\n" "$1" "$n" "$thr" "$nsig"
    sort -t$'\t' -g -k4,4 <(tail -n +2 "pyseer_${tag}.tsv") | awk 'NR<=5' | cut -f1-5
    echo
}
run EPTB  "$SCRIPTS/phenotypes.tsv" 0.05
for d in INH_R RIF_R EMB_R SM_R; do run "$d" "$SCRIPTS/phenotypes.tsv" 0.01; done

echo "===== 5. EPTB: unfiltered vs filtered, side by side"
printf "%-34s %-14s %-14s\n" variant p_unfiltered p_filtered
for v in group_1905 group_1965 group_1946 group_589 group_2137; do
    pu=$(awk -F'\t' -v v="$v" '$1==v{print $4}' "$PY/pyseer_EPTB_genes.tsv" 2>/dev/null | head -1)
    pf=$(awk -F'\t' -v v="$v" '$1 ~ "^"v {print $4}' "pyseer_EPTB_genes_filtered.tsv" 2>/dev/null | head -1)
    printf "%-34s %-14s %-14s\n" "$v" "${pu:-not tested}" "${pf:-removed/not tested}"
done

echo; echo "===== 6. gene-graph neighbourhood of the cas locus and EPTB candidates"
GML=$PAN/final_graph.gml
if [ -s "$GML" ]; then
    for g in csm6 cas1 group_1965 group_1905; do
        echo "--- $g"
        panaroo-gene-neighbourhood --gene "$g" --graph "$GML" --expand_no 4 \
            --out "neighbourhood_${g}.txt" 2>/dev/null \
            && head -12 "neighbourhood_${g}.txt" || echo "    (gene not in graph under this name)"
        echo
    done
else
    echo "  final_graph.gml not found at $GML - skipping"
fi

echo "===== done ($(date)); results: $OUT"