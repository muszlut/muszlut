#!/bin/bash
#SBATCH --job-name=Mtb_7_BioProjects
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/ma95362/Bio_project_publication/Mtb_Multiple_BioProjects/logs/log.%j.out
#SBATCH --error=/scratch/ma95362/Bio_project_publication/Mtb_Multiple_BioProjects/logs/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

set -euo pipefail

# ============================================================
# 1. Activate Bactopia environment
# ============================================================
module load Mamba/23.11.0-0
source activate bactopia4
export PYTHONNOUSERSITE=1

# ============================================================
# 2. Define directories
# ============================================================
BASE_DIR="/scratch/ma95362/Bio_project_publication/Mtb_Multiple_BioProjects"

SEARCH_DIR="${BASE_DIR}/search"
ACCESSION_DIR="${BASE_DIR}/accessions"
RESULT_DIR="${BASE_DIR}/results"
WORK_DIR="${BASE_DIR}/work"

mkdir -p "$SEARCH_DIR" "$ACCESSION_DIR" "$RESULT_DIR" "$WORK_DIR"
cd "$BASE_DIR"

COMBINED="${ACCESSION_DIR}/combined-accessions.txt"

# ============================================================
# 3. BioProject accession list (PRJNA975069 = assemblies only,
#    handled by bactopia_PRJNA975069_assemblies.sh)
# ============================================================
PROJECTS=(
    PRJNA1104194
    PRJNA1204469
    PRJNA935744
    PRJNA1221372
    PRJNA1174701
    PRJNA1182313
    PRJNA1242066
)

echo "=========================================="
echo "Bactopia version: $(bactopia --version)"
echo "=========================================="

# ============================================================
# 4. Search + combine accessions
#    Skipped on resubmission if the combined list already exists,
#    so the input stays identical and -resume can reuse finished tasks.
# ============================================================
if [[ -s "$COMBINED" ]]; then

    echo "Using existing accession list: $COMBINED"

else

    for PROJECT in "${PROJECTS[@]}"; do
        echo "Searching BioProject: ${PROJECT}"
        PROJECT_DIR="${SEARCH_DIR}/${PROJECT}"
        mkdir -p "$PROJECT_DIR"
        cd "$PROJECT_DIR"

        if ! bactopia search --query "$PROJECT" --use-ncbi-genome-size \
             || [[ ! -s bactopia-accessions.txt ]]; then
            echo "WARNING: No accessions for ${PROJECT}; skipping"
            continue
        fi
        echo "Completed: ${PROJECT}"
    done

    cd "$BASE_DIR"
    TEMP="${ACCESSION_DIR}/all-accessions.tmp"
    : > "$TEMP"

    FIRST=1
    for PROJECT in "${PROJECTS[@]}"; do
        F="${SEARCH_DIR}/${PROJECT}/bactopia-accessions.txt"
        [[ -s "$F" ]] || continue
        if [[ $FIRST -eq 1 ]]; then
            head -n 1 "$F" > "$COMBINED"      # header from first file
            FIRST=0
        fi
        tail -n +2 "$F" >> "$TEMP"
    done

    # Remove duplicate Experiment accessions
    awk -F'\t' '!seen[$1]++' "$TEMP" >> "$COMBINED"
    rm -f "$TEMP"
fi

TOTAL=$(awk 'END {print NR-1}' "$COMBINED")
echo "=========================================="
echo "Total unique Experiment accessions: $TOTAL"
echo "=========================================="
[[ "$TOTAL" -gt 0 ]] || { echo "ERROR: No samples retrieved."; exit 1; }

# ============================================================
# 5. Run Bactopia
#    --max_memory below the Slurm --mem (120G) so parallel tasks
#    can never exceed the allocation and get killed by Slurm.
#    Profile kept as 'standard' (conda) because the first 205
#    samples were run with it; switching profile would redo them.
#    -resume skips everything already finished.
# ============================================================
echo "Starting Bactopia analysis..."
cd "$BASE_DIR"

if bactopia \
    --accessions "$COMBINED" \
    --coverage 100 \
    --max_cpus 8 \
    --max_memory 100.GB \
    --max_time 12.h \
    --outdir "$RESULT_DIR" \
    -w "$WORK_DIR" \
    -resume
then
    echo "=========================================="
    echo "Bactopia completed successfully!"
    echo "Total unique samples: $TOTAL"
    echo "Results: $RESULT_DIR"
    echo "=========================================="
else
    echo "=========================================="
    echo "Bactopia failed! Resubmit this script to resume."
    echo "=========================================="
    exit 1
fi