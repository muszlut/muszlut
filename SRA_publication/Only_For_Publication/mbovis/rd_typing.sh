#!/bin/bash
#SBATCH --job-name=rd_typing
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/ma95362/rd_typing/logs/rd.%j.out
#SBATCH --error=/scratch/ma95362/rd_typing/logs/rd.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# =============================================================================
# Region of Difference (RD) typing for the 38 MTBC isolates.
#
# Two independent parts, because the two reference genomes can answer
# different questions:
#
#   PART A  vs H37Rv (TB-Profiler BAMs)
#           RD9, RD4, RD1 are PRESENT in H37Rv and deleted in various MTBC
#           members, so they can only be scored against H37Rv:
#               RD9 deleted, RD4 intact   -> M. caprae
#               RD9 deleted, RD4 deleted  -> M. bovis
#               RD9 intact                -> M. tuberculosis sensu stricto
#               RD1 also deleted          -> BCG
#           Coordinates are looked up from the H37Rv GFF by Rv locus tag,
#           not hardcoded, so they are correct for this assembly.
#
#   PART B  vs M. bovis AF2122/97 (vSNP3 BAMs)
#           RDAf2 is deleted in African 2 strains and PRESENT in AF2122/97
#           (a Eu1 strain), so it is scored here. Rather than trusting a
#           published coordinate, this does a data-driven scan: it finds
#           windows with no coverage in the 20 Af2CC samples that ARE
#           covered in the 14 non-Af2 M. bovis samples. A region that
#           partitions the dataset exactly along the Af2 split is RDAf2.
#
# Run: sbatch rd_typing.sh
# =============================================================================

module load SAMtools

set -euo pipefail

# ----------------------------- configuration ---------------------------------

OUT="/scratch/ma95362/rd_typing"
BACTOPIA="/scratch/ma95362/Bio_project_publication/Mbovis"
VSNP="/scratch/ma95362/vsnp3/step1"

H37RV_GFF="/scratch/ma95362/gbk/ncbi_dataset/data/GCF_000195955.2/genomic.gff"
WINDOW=1000          # bp, for the part B deletion scan
MIN_DEPTH=5          # a window below this mean depth counts as absent

# Group membership from the vSNP3 step 2 output
AF2="17 151P 213P 262G 264G 264P 271G 281G 281P 283G 292G 300P 301G 304G 304P 305G 308P 310P 342P 343P"
NONAF2="120PG 231G 258G 259G 262P 288P 307G 318G 322G 323G 32_HK 341P 147 234GP"
CAPRAE="097G_400 100GP 106GP 97G"

mkdir -p "$OUT/logs"
cd "$OUT"

echo "Started: $(date)"
echo

# =============================================================================
# PART A - RD9 / RD4 / RD1 against H37Rv
# =============================================================================

echo "############ PART A: RD typing against H37Rv ############"
echo

if [ ! -s "$H37RV_GFF" ]; then
    echo "SKIPPING PART A: H37Rv GFF not found at $H37RV_GFF" >&2
    echo "Download with:" >&2
    echo "  datasets download genome accession GCF_000195955.2 --include genome,gff3,gbff" >&2
else

# Look up the flanking genes of each RD by Rv locus tag and take the span.
# RD9  : Rv2072c - Rv2075c
# RD4  : Rv1505c - Rv1516c
# RD1  : Rv3871  - Rv3879c
span_from_tags () {
    local first="$1" last="$2"
    awk -F'\t' -v a="$first" -v b="$last" '
        $3=="gene" && $9 ~ ("locus_tag=" a "(;|$)") { s=$4 }
        $3=="gene" && $9 ~ ("locus_tag=" b "(;|$)") { e=$5 }
        END { if (s && e) { if (s>e) {t=s; s=e; e=t}; print s"\t"e } }
    ' "$H37RV_GFF"
}

PROBE=$(ls "$BACTOPIA"/*/tools/tbprofiler/supplemental/*.bam 2>/dev/null | head -1)
CHROM=$(samtools view -H "$PROBE" | awk '$1=="@SQ"{sub("SN:","",$2); print $2; exit}')
echo "H37Rv contig (from BAM): $CHROM"

declare -A RD
for spec in "RD9:Rv2072c:Rv2075c" "RD4:Rv1505c:Rv1516c" "RD1:Rv3871:Rv3879c"; do
    name=${spec%%:*}; rest=${spec#*:}
    first=${rest%%:*}; last=${rest##*:}
    coords=$(span_from_tags "$first" "$last")
    if [ -n "$coords" ]; then
        RD[$name]="$coords"
        printf "  %-5s %s..%s  (%s - %s)\n" "$name" $(echo "$coords") "$first" "$last"
    else
        echo "  $name: could not resolve $first / $last in the GFF" >&2
    fi
done
echo

{
printf "sample"
for n in RD9 RD4 RD1; do printf "\t%s_mean\t%s_pct_zero" "$n" "$n"; done
printf "\tcall\n"

for s in $AF2 $NONAF2 $CAPRAE; do
    b="$BACTOPIA/$s/tools/tbprofiler/supplemental/$s.bam"
    [ -f "$b" ] || { echo "$s: no H37Rv bam" >&2; continue; }
    [ -f "$b.bai" ] || samtools index "$b" 2>/dev/null || true

    printf "%s" "$s"
    declare -A pz
    for n in RD9 RD4 RD1; do
        [ -n "${RD[$n]:-}" ] || { printf "\tNA\tNA"; pz[$n]=-1; continue; }
        read -r st en <<< "${RD[$n]}"
        read -r mean pct <<< "$(samtools depth -a -r "${CHROM}:${st}-${en}" "$b" \
            | awk '{t+=$3; n++; if($3==0) z++} END {
                     if(n==0) print "0 100"; else printf "%.1f %.1f", t/n, 100*(z+0)/n }')"
        printf "\t%s\t%s" "$mean" "$pct"
        pz[$n]=${pct%%.*}
    done

    # interpret: >80% of the region at zero coverage = deleted
    call="unresolved"
    if   [ "${pz[RD9]}" -gt 80 ] && [ "${pz[RD4]}" -gt 80 ]; then call="M. bovis"
    elif [ "${pz[RD9]}" -gt 80 ] && [ "${pz[RD4]}" -le 80 ]; then call="M. caprae"
    elif [ "${pz[RD9]}" -le 80 ];                            then call="M. tuberculosis s.s."
    fi
    [ "${pz[RD1]}" -gt 80 ] && call="$call + RD1-deleted (BCG-like)"

    printf "\t%s\n" "$call"
done
} | tee "$OUT/rd_typing_h37rv.tsv" | column -t -s$'\t'

fi

echo
echo "Part A written to $OUT/rd_typing_h37rv.tsv"
echo

# =============================================================================
# PART B - data-driven RDAf2 scan against AF2122/97
# =============================================================================

echo "############ PART B: deletion scan against AF2122/97 ############"
echo "Looking for windows absent in Af2CC and present in non-Af2 M. bovis."
echo

mkdir -p "$OUT/cov"

# Mean depth per window, per sample
for s in $AF2 $NONAF2; do
    cf="$OUT/cov/${s}.win.txt"
    [ -s "$cf" ] && continue
    b="$VSNP/$s/alignment_NC_002945v4/${s}_nodup.bam"
    [ -f "$b" ] || { echo "$s: no AF2122 bam" >&2; continue; }
    echo "  scanning $s"
    samtools depth -a "$b" \
      | awk -v w="$WINDOW" '{ bin=int($2/w); t[bin]+=$3; n[bin]++ }
             END { for (b in t) printf "%d\t%.1f\n", b*w, t[b]/n[b] }' \
      | sort -k1,1n > "$cf"
done

echo
echo "  comparing groups..."

python3 - "$OUT" "$WINDOW" "$MIN_DEPTH" "$AF2" "$NONAF2" <<'PY'
import sys, os, glob
out, win, mind = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
af2    = sys.argv[4].split()
nonaf2 = sys.argv[5].split()

def load(s):
    p = os.path.join(out, "cov", f"{s}.win.txt")
    if not os.path.exists(p): return None
    d = {}
    for line in open(p):
        a, b = line.split()
        d[int(a)] = float(b)
    return d

cov = {s: load(s) for s in af2 + nonaf2}
cov = {s: v for s, v in cov.items() if v}
A = [s for s in af2 if s in cov]
B = [s for s in nonaf2 if s in cov]
print(f"  Af2 samples: {len(A)}   non-Af2 samples: {len(B)}")

positions = set()
for v in cov.values(): positions |= set(v)

hits = []
for p in sorted(positions):
    a = [cov[s].get(p, 0) for s in A]
    b = [cov[s].get(p, 0) for s in B]
    if not a or not b: continue
    if max(a) < mind and min(b) >= mind:
        hits.append(p)

# merge adjacent windows into regions
regions, start, prev = [], None, None
for p in hits:
    if start is None:
        start = prev = p
    elif p == prev + win:
        prev = p
    else:
        regions.append((start, prev + win)); start = prev = p
if start is not None: regions.append((start, prev + win))

with open(os.path.join(out, "af2_specific_deletions.tsv"), "w") as fh:
    fh.write("start\tend\tlength_bp\n")
    print("\n  Regions absent in ALL Af2 and present in ALL non-Af2:")
    if not regions:
        print("    (none found)")
    for s, e in regions:
        fh.write(f"{s}\t{e}\t{e-s}\n")
        print(f"    NC_002945.4:{s}-{e}   {e-s} bp")
PY

echo
echo "Part B written to $OUT/af2_specific_deletions.tsv"
echo
echo "Annotate any region found with:"
echo "  awk -F'\\t' -v s=<start> -v e=<end> '\$3==\"gene\" && \$4<e && \$5>s' \\"
echo "    /scratch/ma95362/ref_mbovis/ncbi_dataset/data/GCF_000195835.3/genomic.gff"
echo
echo "A single region of roughly 10-15 kb that splits the dataset exactly along"
echo "the Af2 / non-Af2 line is RDAf2. Check its genes against Berg et al. 2011"
echo "(J Bacteriol, African 2 clonal complex) before naming it in the paper."
echo
echo "Finished: $(date)"