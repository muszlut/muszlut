#!/bin/bash
#SBATCH --job-name=fastq_recover
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/ma95362/Bio_project_publication/Mbovis/logs/recover.%j.out
#SBATCH --error=/scratch/ma95362/Bio_project_publication/Mbovis/logs/recover.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# =============================================================================
# Recover readable reads from CRC-corrupted gzipped FASTQ files.
#
#   Step 1  salvage  - decompress each corrupt .fastq.gz until gzip hits the
#                      corruption point, keeping only COMPLETE 4-line records.
#   Step 2  pair     - R1 and R2 break at different points, so restore sync.
#   Step 3  verify   - read counts, pairing, gzip integrity, depth estimate.
#
# The original files are never modified.
#
# ---------------------------------------------------------------------------
# RUN HISTORY
#   job 49892221  failed at preflight (seqkit not loaded) - no work done
#   job 49892389  step 1 OK, step 2 killed at exactly 16.0 GB (MaxRSS 16775636K)
#                 seqkit pair loads the whole first file into memory; the kill
#                 left truncated gzip streams ("unexpected end of file").
#
# Step 1 output is verified clean, so STEP 1 IS NOW DISABLED (RUN_SALVAGE=0).
# Memory raised 16G -> 64G for the re-run of step 2.
#
# To re-run step 1 from scratch (e.g. fresh source files), set RUN_SALVAGE=1.
# =============================================================================

RUN_SALVAGE=0             # 1 = run step 1 salvage, 0 = skip (already done)
RUN_PAIR=1                # 1 = run step 2 pairing
RUN_VERIFY=1              # 1 = run step 3 verification

# --------------------------- environment -------------------------------------
# A batch shell does not inherit the interactive conda environment.

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

if command -v pigz >/dev/null 2>&1; then
    ZIP="pigz -p ${SLURM_CPUS_PER_TASK:-4}"
else
    ZIP="gzip"
fi

[ -d "$SRC" ] || { echo "ERROR: source directory not found: $SRC" >&2; exit 1; }
mkdir -p "$OUT"

# Free space guard: salvage + paired output need roughly 50 GB together.
avail_kb=$(df -P "$OUT" | awk 'NR==2 {print $4}')
if [ "$avail_kb" -lt 60000000 ]; then
    echo "ERROR: only $((avail_kb/1024/1024)) GB free on $OUT; need ~60 GB." >&2
    exit 1
fi

REPORT="$OUT/recovery_report.tsv"

echo "Source:      $SRC"
echo "Destination: $OUT"
echo "Compressor:  $ZIP"
echo "Memory:      ${SLURM_MEM_PER_NODE:-unknown} MB requested"
echo "Steps:       salvage=$RUN_SALVAGE pair=$RUN_PAIR verify=$RUN_VERIFY"
echo "Started:     $(date)"
echo

# ========================== STEP 1: SALVAGE PREFIX ===========================
# DISABLED by default - job 49892389 already produced verified-clean output.
# Single pass per file: decompress, emit only whole 4-line records, recompress.

if [ "$RUN_SALVAGE" -eq 1 ]; then

    printf "sample\tread\tcomplete_records\tbytes_readable\n" > "$REPORT"

    for s in $BAD; do
        for r in R1 R2; do
            f="${s}_${r}.fastq.gz"

            if [ ! -f "$SRC/$f" ]; then
                echo "WARNING: missing $f -- skipping" >&2
                continue
            fi

            echo "[salvage] $f"

            # awk buffers four lines and only emits a record once all four
            # exist, so a partial record at the corruption point is dropped.
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

else
    echo "[step 1] SKIPPED (RUN_SALVAGE=0) - using existing salvaged files."
    echo
fi

# ============================ STEP 2: RESTORE PAIRING ========================
# R1 and R2 stopped at different records, so the two files are out of sync.
# seqkit pair keeps only reads present in both.
#
# Changes after the job 49892389 failure:
#   - memory raised to 64G (seqkit holds the whole first file in RAM)
#   - salvaged inputs are gzip-tested BEFORE pairing, so a bad input aborts
#     the run instead of silently producing garbage
#   - seqkit's exit status is now checked; a failure stops the script
#   - stale truncated outputs from the previous run are removed first

if [ "$RUN_PAIR" -eq 1 ]; then

    cd "$OUT" || exit 1

    # -- 2a. validate the salvaged inputs before spending time on them --------
    echo "[check] verifying salvaged inputs before pairing"
    bad_input=0
    for s in $BAD; do
        for r in R1 R2; do
            f="$OUT/${s}_${r}.fastq.gz"
            if [ ! -f "$f" ]; then
                echo "  MISSING: ${s}_${r}.fastq.gz" >&2
                bad_input=1
            elif ! gzip -t "$f" 2>/dev/null; then
                echo "  TRUNCATED: ${s}_${r}.fastq.gz" >&2
                bad_input=1
            fi
        done
    done
    if [ "$bad_input" -eq 1 ]; then
        echo "ERROR: salvaged inputs are not all valid. Set RUN_SALVAGE=1" >&2
        echo "       and re-run to regenerate them." >&2
        exit 1
    fi
    echo "        all $(echo $BAD | wc -w) samples OK"
    echo

    # -- 2b. clear stale output from the killed run --------------------------
    mkdir -p "$OUT/paired"
    rm -f "$OUT/paired"/*.fastq.gz

    # -- 2c. pair ------------------------------------------------------------
    for s in $BAD; do
        if [ ! -f "${s}_R1.fastq.gz" ] || [ ! -f "${s}_R2.fastq.gz" ]; then
            echo "WARNING: $s is missing a mate file -- skipping pairing" >&2
            continue
        fi

        echo "[pair] $s"
        if ! seqkit pair -1 "${s}_R1.fastq.gz" -2 "${s}_R2.fastq.gz" \
                         -u -O paired 2>&1 | sed 's/^/           /'; then
            echo "ERROR: seqkit pair failed for $s." >&2
            echo "       If this is an out-of-memory kill, raise --mem further" >&2
            echo "       or switch to the awk fallback described in the notes." >&2
            exit 1
        fi

        # confirm this sample's output is a complete gzip stream before moving on
        for r in R1 R2; do
            if ! gzip -t "$OUT/paired/${s}_${r}.fastq.gz" 2>/dev/null; then
                echo "ERROR: ${s}_${r} paired output is truncated." >&2
                echo "       Almost certainly an out-of-memory kill. Raise --mem." >&2
                exit 1
            fi
        done
    done

    echo
    echo "Pairing complete: $(date)"
    echo

else
    echo "[step 2] SKIPPED (RUN_PAIR=0)"
    echo
fi

# ============================== STEP 3: VERIFY ===============================

if [ "$RUN_VERIFY" -eq 1 ]; then

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

        if gzip -t "$f1" 2>/dev/null && gzip -t "$f2" 2>/dev/null; then
            gzip_ok="OK"
        else
            gzip_ok="BAD"
        fi

        a=$(( $(zcat "$f1" 2>/dev/null | wc -l) / 4 ))
        b=$(( $(zcat "$f2" 2>/dev/null | wc -l) / 4 ))

        [ "$a" -eq "$b" ] && pair_ok="OK" || pair_ok="MISMATCH"

        depth=$(( (a + b) * READ_LEN / GENOME_SIZE ))

        if [ "$gzip_ok" != "OK" ]; then
            verdict="FAIL_TRUNCATED"
        elif [ "$pair_ok" != "OK" ]; then
            verdict="FAIL_UNPAIRED"
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

    usable=$(awk -F'\t' 'NR>1 && $7=="USABLE"' "$SUMMARY" | wc -l)
    total=$(echo $BAD | wc -w)
    echo "$usable of $total samples are USABLE."
    echo

    echo "Per-file salvage log:  $REPORT"
    echo "Verification summary:  $SUMMARY"
    echo "Recovered FASTQs:      $OUT/paired/"
    echo
    echo "Depth is estimated as (R1+R2 reads x ${READ_LEN} bp) / ${GENOME_SIZE} bp."
    echo "Confirm your real read length before trusting it:"
    echo "  zcat $OUT/paired/343P_R1.fastq.gz | head -2 | tail -1 | wc -c"
    echo
    echo "Next: feed $OUT/paired/ to 'bactopia prepare' for these samples, and"
    echo "      keep the 49 clean samples pointing at $SRC."
fi

echo
echo "Finished: $(date)"