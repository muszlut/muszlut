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
# All coordinates below were determined EMPIRICALLY from these BAMs rather
# than taken from the literature, by finding the exact zero-coverage
# boundaries and confirming they are identical across every member of a
# group. Verify the gene content against the published RD definitions
# before using the published RD names in the manuscript.
#
# ---------------------------------------------------------------------------
# WHAT WAS FOUND (H37Rv / NC_000962.3, BAM contig name "Chromosome")
#
#   RD9 core        2,330,074 - 2,332,101   (2,028 bp)
#       Deleted in ALL 38 isolates -> none is M. tuberculosis sensu stricto.
#       NB the Rv2072c-Rv2075c span used earlier includes the intact
#       FLANKING genes, which is why every sample read ~51% zero and a
#       >80% threshold never fired.
#
#   M. bovis deletion   1,696,017 - 1,708,748  (12,732 bp)
#       Identical outer boundaries in all 34 M. bovis. Some samples show
#       ~20 bp islands of coverage inside it (e.g. 1,701,629-1,701,647);
#       those are mismapped reads at 1-2x against ~40-90x flanking.
#
#   La2 deletion        1,688,917 - 1,700,496  (11,580 bp)
#       Identical in all four La2 samples (97G, 100GP, 106GP, 097G_400),
#       breakpoints agreeing to within 1 bp.
#
#   The two deletions OVERLAP only between 1,700,497 and 1,696,017, so a
#   simple "RD4 percent deleted" is ambiguous. This script instead scores
#   the two NON-overlapping windows, which give a clean binary:
#       BOV_ONLY  1,700,497 - 1,708,748  deleted in M. bovis, present in La2
#       CAP_ONLY  1,688,917 - 1,696,016  deleted in La2, present in M. bovis
#
#   RD1             4,348,827 - 4,359,782
#       0% zero coverage in all 38 -> RD1 intact, no BCG.
#
# PART B re-runs the data-driven Af2 deletion scan against AF2122/97.
#   Previously found: NC_002945.4:681,580-695,673 (14,094 bp), identical
#   breakpoints in all 20 Af2 isolates, spanning BQ2027_MB0599-MB0610
#   including mce2R. Candidate RDAf2 - confirm against Berg et al. 2011.
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
MIN_DEPTH=5          # window mean below this counts as absent (part B)
DEPTH_CUT=3          # per-base depth below this counts as absent (part A)
ZERO_CALL=80         # % of a window below DEPTH_CUT to call it deleted

AF2="17 151P 213P 262G 264G 264P 271G 281G 281P 283G 292G 300P 301G 304G 304P 305G 308P 310P 342P 343P"
NONAF2="120PG 231G 258G 259G 262P 288P 307G 318G 322G 323G 32_HK 341P 147 234GP"
CAPRAE="097G_400 100GP 106GP 97G"

# Empirically determined intervals (see header)
RD9_CORE="2330074-2332101"
BOV_ONLY="1700497-1708748"
CAP_ONLY="1688917-1696016"
RD1_SPAN="4348827-4359782"

mkdir -p "$OUT/logs" "$OUT/cov"
cd "$OUT"

echo "Started: $(date)"
echo

# =============================================================================
# PART A - RD typing against H37Rv
# =============================================================================

echo "############ PART A: RD typing against H37Rv ############"
echo

PROBE=$(ls "$BACTOPIA"/*/tools/tbprofiler/supplemental/*.bam 2>/dev/null | head -1)
[ -n "$PROBE" ] || { echo "ERROR: no TB-Profiler BAMs found under $BACTOPIA" >&2; exit 1; }
CHROM=$(samtools view -H "$PROBE" | awk '$1=="@SQ"{sub("SN:","",$2); print $2; exit}')
echo "Contig from BAM header: $CHROM"
echo "  RD9 core  $RD9_CORE"
echo "  BOV_ONLY  $BOV_ONLY   (deleted in M. bovis, present in La2)"
echo "  CAP_ONLY  $CAP_ONLY   (deleted in La2, present in M. bovis)"
echo "  RD1       $RD1_SPAN"
echo

# mean depth and % of bases below DEPTH_CUT for one region
region_stats () {
    local bam="$1" reg="$2"
    samtools depth -a -r "${CHROM}:${reg}" "$bam" 2>/dev/null \
      | awk -v c="$DEPTH_CUT" '{t+=$3; n++; if($3<c) z++}
             END { if(n==0) print "0 100"; else printf "%.1f %.1f", t/n, 100*(z+0)/n }'
}

{
printf "sample\tgroup\tRD9_mean\tRD9_pctdel\tBOVonly_mean\tBOVonly_pctdel"
printf "\tCAPonly_mean\tCAPonly_pctdel\tRD1_mean\tRD1_pctdel\tcall\n"

for s in $AF2 $NONAF2 $CAPRAE; do
    case " $CAPRAE " in *" $s "*) grp="La2" ;; *) case " $AF2 " in *" $s "*) grp="Af2" ;; *) grp="nonAf2" ;; esac ;; esac

    b="$BACTOPIA/$s/tools/tbprofiler/supplemental/$s.bam"
    [ -f "$b" ] || { echo "$s: no H37Rv bam" >&2; continue; }

    read -r rd9_m  rd9_z  <<< "$(region_stats "$b" "$RD9_CORE")"
    read -r bov_m  bov_z  <<< "$(region_stats "$b" "$BOV_ONLY")"
    read -r cap_m  cap_z  <<< "$(region_stats "$b" "$CAP_ONLY")"
    read -r rd1_m  rd1_z  <<< "$(region_stats "$b" "$RD1_SPAN")"

    rd9i=${rd9_z%%.*}; bovi=${bov_z%%.*}; capi=${cap_z%%.*}; rd1i=${rd1_z%%.*}

    if [ "$rd9i" -lt "$ZERO_CALL" ]; then
        call="M. tuberculosis sensu stricto (RD9 intact)"
    elif [ "$bovi" -ge "$ZERO_CALL" ] && [ "$capi" -lt "$ZERO_CALL" ]; then
        call="M. bovis"
    elif [ "$capi" -ge "$ZERO_CALL" ] && [ "$bovi" -lt "$ZERO_CALL" ]; then
        call="M. caprae-like (La2)"
    else
        call="unresolved"
    fi
    [ "$rd1i" -ge "$ZERO_CALL" ] && call="$call; RD1 DELETED (BCG-like)"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$s" "$grp" "$rd9_m" "$rd9_z" "$bov_m" "$bov_z" \
        "$cap_m" "$cap_z" "$rd1_m" "$rd1_z" "$call"
done
} | tee "$OUT/rd_typing_h37rv.tsv" | column -t -s$'\t'

echo
echo "Part A -> $OUT/rd_typing_h37rv.tsv"
echo

# ------------------------- annotate the two intervals ------------------------

if [ -s "$H37RV_GFF" ]; then
  for spec in "M.bovis_deletion:1696017:1708748" "La2_deletion:1688917:1700496" "RD9_core:2330074:2332101"; do
      n=${spec%%:*}; r=${spec#*:}; st=${r%%:*}; en=${r##*:}
      echo "=== genes in $n ($st-$en)"
      awk -F'\t' -v s="$st" -v e="$en" '$3=="gene" && $4<e && $5>s {
          match($0,/locus_tag=[^;]*/); lt=substr($0,RSTART+10,RLENGTH-10)
          nm="-"; if (match($0,/;Name=[^;]*/)) nm=substr($0,RSTART+6,RLENGTH-6)
          printf "   %8d-%-8d %-16s %s\n", $4, $5, lt, nm }' "$H37RV_GFF"
      echo
  done | tee "$OUT/rd_gene_content.txt"
fi

# =============================================================================
# PART B - data-driven Af2 deletion scan against AF2122/97
# =============================================================================

echo "############ PART B: deletion scan against AF2122/97 ############"
echo

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

python3 - "$OUT" "$WINDOW" "$MIN_DEPTH" "$AF2" "$NONAF2" <<'PY'
import sys, os
out, win, mind = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
af2, nonaf2 = sys.argv[4].split(), sys.argv[5].split()

def load(s):
    p = os.path.join(out, "cov", f"{s}.win.txt")
    if not os.path.exists(p): return None
    return {int(a): float(b) for a, b in (l.split() for l in open(p))}

cov = {s: v for s, v in ((s, load(s)) for s in af2 + nonaf2) if v}
A = [s for s in af2 if s in cov]; B = [s for s in nonaf2 if s in cov]
print(f"  Af2: {len(A)}   non-Af2: {len(B)}")

positions = set().union(*(set(v) for v in cov.values()))
hits = [p for p in sorted(positions)
        if max(cov[s].get(p, 0) for s in A) < mind
        and min(cov[s].get(p, 0) for s in B) >= mind]

regions, start, prev = [], None, None
for p in hits:
    if start is None: start = prev = p
    elif p == prev + win: prev = p
    else: regions.append((start, prev + win)); start = prev = p
if start is not None: regions.append((start, prev + win))

with open(os.path.join(out, "af2_specific_deletions.tsv"), "w") as fh:
    fh.write("start\tend\tlength_bp\n")
    print("\n  Absent in ALL Af2, present in ALL non-Af2:")
    for s, e in regions or []:
        fh.write(f"{s}\t{e}\t{e-s}\n")
        print(f"    NC_002945.4:{s}-{e}   {e-s} bp")
    if not regions: print("    (none)")
PY

echo
echo "Part B -> $OUT/af2_specific_deletions.tsv"
echo
echo "Previously refined to base resolution: NC_002945.4:681,580-695,673"
echo "(14,094 bp), identical in all 20 Af2 isolates, spanning"
echo "BQ2027_MB0599-MB0610 including mce2R. Candidate RDAf2 -- confirm the"
echo "gene content against Berg et al. 2011, J Bacteriol 193(3):670."
echo
echo "Finished: $(date)"