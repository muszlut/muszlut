#!/bin/bash
#SBATCH --job-name=Mtb_7_BioProjects_sing
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

# Containers, same as the Eth763 and PRJNA975069 runs
module load Apptainer 2>/dev/null || module load Singularity 2>/dev/null || true
export NXF_SINGULARITY_CACHEDIR=/scratch/ma95362/singularity_cache
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

# ============================================================
# 2. Define directories
# ============================================================
BASE_DIR="/scratch/ma95362/Bio_project_publication/Mtb_Multiple_BioProjects"

SEARCH_DIR="${BASE_DIR}/search"
ACCESSION_DIR="${BASE_DIR}/accessions"
RESULT_DIR="${BASE_DIR}/results_singularity"   # fresh, separate from the old conda run
WORK_DIR="${BASE_DIR}/work_singularity"

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
#    -profile singularity: the conda ('standard') run kept hanging in QC,
#    while the singularity runs (Eth763, PRJNA975069) worked. Using the
#    same profile for all three datasets also keeps the methods identical.
#    -resume skips everything already finished if you resubmit.
# ============================================================
echo "Starting Bactopia analysis..."
cd "$BASE_DIR"

# ------------------------------------------------------------
# Custom Nextflow config:
#  - time: Bactopia's QC step requests 2 h by default and
#    --max_time only caps that, so set the time on every step directly.
#  - errorStrategy: retry a failed sample once, then skip it so
#    one bad sample does not kill the whole run.
#    Skipped samples are listed in results/bactopia-runs/*/nf-reports.
# Changing these settings does not invalidate -resume.
# ------------------------------------------------------------
CUSTOM_CONFIG="${BASE_DIR}/mtb_custom.config"
cat > "$CUSTOM_CONFIG" <<'EOF'
process {
    withName: '.*' {
        time          = '12.h'
        errorStrategy = { task.attempt <= 1 ? 'retry' : 'ignore' }
        maxRetries    = 1
    }
}
EOF

if bactopia \
    --accessions "$COMBINED" \
    --coverage 100 \
    --max_cpus 8 \
    --max_memory 100.GB \
    --max_time 12.h \
    --outdir "$RESULT_DIR" \
    -w "$WORK_DIR" \
    -c "$CUSTOM_CONFIG" \
    -profile singularity \
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