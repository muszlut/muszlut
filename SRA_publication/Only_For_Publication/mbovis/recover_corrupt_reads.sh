#!/bin/bash
#SBATCH --job-name=fastq_recover
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/ma95362/Bio_project_publication/Mbovis/logs/recover.%j.out
#SBATCH --error=/scratch/ma95362/Bio_project_publication/Mbovis/logs/recover.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# =============================================================================
# Recover readable reads from CRC-corrupted gzipped FASTQ files.
#
# Each corrupt .fastq.gz is decompressed until gzip hits the corruption point.
# Only COMPLETE 4-line records from the readable prefix are kept. R1 and R2
# break at different points, so the files are then re-paired to restore sync.
#
# The original files are never modified.
#
# Run:  sbatch recover_corrupt_reads.sh
# =============================================================================

# --------------------------- environment -------------------------------------
# A batch shell does not inherit the interactive conda environment, so load
# what we need here. SeqKit comes from the cluster module system.

module load SeqKit/2.9.0

set -uo pipefail          # NOTE: deliberately no -e; zcat exits non-zero on
                          # every corrupt file and that is expected here.

# ----------------------------- configuration ---------------------------------

SRC="/scratch/ma95362/Bio_project_publicatio/mbovis_seq"
OUT="/scratch/ma95362/Bio_project_publicatio/mbovis_recovered"
READ_LEN=150              # only used for the depth estimate in the report
GENOME_SIZE=4345492       # M. bovis AF2122/97
MIN_DEPTH=30              # samples below this are flagged as unusable

BAD="288P 305G 318G 320P 322G 323G 32_HK 336G 341P 342P 343P 379 380 383 SU109"

# ------------------------------ preflight ------------------------------------

if ! command -v seqkit >/dev/null 2>&1; then
    echo "ERROR: seqkit not on PATH after 'module load SeqKit/2.9.0'." >&2
    echo "       Check 'module spider SeqKit' for an available version." >&2
    exit 1
fi
echo "seqkit: $(command -v seqkit) ($(seqkit version 2>&1 | head -1))"

# pigz gives a large speedup on recompression if it is available
if command -v pigz >/dev/null 2>&1; then
    ZIP="pigz -p ${SLURM_CPUS_PER_TASK:-4}"
else
    ZIP="gzip"
fi

[ -d "$SRC" ] || { echo "ERROR: source directory not found: $SRC" >&2; exit 1; }
mkdir -p "$OUT"

REPORT="$OUT/recovery_report.tsv"
printf "sample\tread\tcomplete_records\tbytes_readable\n" > "$REPORT"

echo "Source:      $SRC"
echo "Destination: $OUT"
echo "Compressor:  $ZIP"
echo "Started:     $(date)"
echo

# -------------------------- step 1: salvage prefix ---------------------------
# Single pass per file: decompress, emit only whole 4-line records, recompress.
# Doing this in one pass rather than counting lines first halves the I/O.

for s in $BAD; do
    for r in R1 R2; do
        f="${s}_${r}.fastq.gz"

        if [ ! -f "$SRC/$f" ]; then
            echo "WARNING: missing $f -- skipping" >&2
            continue
        fi

        echo "[salvage] $f"

        # awk buffers four lines and only emits a record once all four exist,
        # so a partial record at the corruption point is silently dropped.
        # The record count is written to a side file rather than stderr so it
        # survives the pipeline cleanly.
        countfile="$OUT/.count.$$"
        zcat "$SRC/$f" 2>/dev/null \
          | awk -v cf="$countfile" '
                { buf[NR % 4] = $0 }
                NR % 4 == 0 {
                    print buf[1]; print buf[2]; print buf[3]; print buf[0]
                    kept++
                }
                END { print kept + 0 > cf }
            ' \
          | $ZIP > "$OUT/$f"

        n=$(cat "$countfile" 2>/dev/null || echo 0)
        rm -f "$countfile"

        bytes=$(stat -c %s "$OUT/$f" 2>/dev/null || echo 0)
        printf "%s\t%s\t%s\t%s\n" "$s" "$r" "${n:-0}" "$bytes" >> "$REPORT"
        echo "           recovered ${n:-0} records"
    done
done

echo
echo "Salvage complete: $(date)"
echo

# --------------------------- step 2: restore pairing -------------------------
# R1 and R2 almost certainly stopped at different records, so the two files
# are out of sync. seqkit pair keeps only reads present in both and writes
# them to the directory given by -O.

cd "$OUT" || exit 1

for s in $BAD; do
    if [ ! -f "${s}_R1.fastq.gz" ] || [ ! -f "${s}_R2.fastq.gz" ]; then
        echo "WARNING: $s is missing a mate file -- skipping pairing" >&2
        continue
    fi
    echo "[pair] $s"
    seqkit pair -1 "${s}_R1.fastq.gz" -2 "${s}_R2.fastq.gz" -u -O paired 2>&1 \
        | sed 's/^/           /'
done

echo
echo "Pairing complete: $(date)"
echo

# ------------------------------ step 3: verify -------------------------------

mkdir -p "$OUT/paired"
SUMMARY="$OUT/paired/verification.tsv"
printf "sample\tR1_reads\tR2_reads\tpairing\tgzip\test_depth_x\tverdict\n" > "$SUMMARY"

echo "=============================================================="
printf "%-8s %12s %12s %8s %8s %9s  %s\n" \
       SAMPLE R1 R2 PAIRED GZIP DEPTH VERDICT
echo "=============================================================="

for s in $BAD; do
    f1="$OUT/paired/${s}_R1.fastq.gz"
    f2="$OUT/paired/${s}_R2.fastq.gz"

    if [ ! -f "$f1" ] || [ ! -f "$f2" ]; then
        printf "%-8s %12s %12s %8s %8s %9s  %s\n" "$s" "-" "-" "-" "-" "-" "MISSING"
        printf "%s\t0\t0\tNA\tNA\t0\tMISSING\n" "$s" >> "$SUMMARY"
        continue
    fi

    a=$(( $(zcat "$f1" | wc -l) / 4 ))
    b=$(( $(zcat "$f2" | wc -l) / 4 ))

    [ "$a" -eq "$b" ] && pair_ok="OK" || pair_ok="MISMATCH"

    if gzip -t "$f1" 2>/dev/null && gzip -t "$f2" 2>/dev/null; then
        gzip_ok="OK"
    else
        gzip_ok="BAD"
    fi

    depth=$(( (a + b) * READ_LEN / GENOME_SIZE ))

    if [ "$pair_ok" != "OK" ] || [ "$gzip_ok" != "OK" ]; then
        verdict="FAIL"
    elif [ "$depth" -lt "$MIN_DEPTH" ]; then
        verdict="LOW_DEPTH"
    else
        verdict="USABLE"
    fi

    printf "%-8s %12s %12s %8s %8s %8sx  %s\n" \
           "$s" "$a" "$b" "$pair_ok" "$gzip_ok" "$depth" "$verdict"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
           "$s" "$a" "$b" "$pair_ok" "$gzip_ok" "$depth" "$verdict" >> "$SUMMARY"
done

echo "=============================================================="
echo
echo "Per-file salvage log:  $REPORT"
echo "Verification summary:  $SUMMARY"
echo "Recovered FASTQs:      $OUT/paired/"
echo
echo "Depth is estimated as (R1+R2 reads x ${READ_LEN} bp) / ${GENOME_SIZE} bp."
echo "Confirm your real read length before trusting it:"
echo "  zcat $OUT/paired/343P_R1.fastq.gz | head -2 | tail -1 | wc -c"
echo
echo "Next: feed $OUT/paired/ to 'bactopia prepare' for these 15 samples, and"
echo "      keep the 49 clean samples pointing at $SRC."
echo
echo "Finished: $(date)"