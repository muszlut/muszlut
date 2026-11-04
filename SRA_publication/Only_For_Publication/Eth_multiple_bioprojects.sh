#!/bin/bash
#SBATCH --job-name=Mtb_8_BioProjects
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

# Load Mamba
module load Mamba/23.11.0-0

# Activate Bactopia environment
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

mkdir -p "$SEARCH_DIR"
mkdir -p "$ACCESSION_DIR"
mkdir -p "$RESULT_DIR"
mkdir -p "$WORK_DIR"

cd "$BASE_DIR"

# ============================================================
# 3. BioProject accession list
# ============================================================

PROJECTS=(
    PRJNA1104194
    PRJNA1204469
    PRJNA975059
    PRJNA935744
    PRJNA1221372
    PRJNA1174701
    PRJNA1182212
    PRJNA1242066
)

# ============================================================
# 4. Check Bactopia version
# ============================================================

echo "=========================================="
echo "Bactopia version:"
bactopia --version
echo "=========================================="

# ============================================================
# 5. Search all eight BioProjects
# ============================================================

for PROJECT in "${PROJECTS[@]}"; do

    echo "Searching BioProject: ${PROJECT}"

    PROJECT_DIR="${SEARCH_DIR}/${PROJECT}"

    mkdir -p "$PROJECT_DIR"

    cd "$PROJECT_DIR"

    bactopia search \
        --query "$PROJECT" \
        --use-ncbi-genome-size

    if [[ ! -s bactopia-accessions.txt ]]; then

        echo "ERROR: No accessions found for ${PROJECT}"

        exit 1

    fi

    echo "Completed: ${PROJECT}"

done

# ============================================================
# 6. Combine all Experiment accessions
# ============================================================

echo "Combining all BioProject accessions..."

cd "$BASE_DIR"

COMBINED="${ACCESSION_DIR}/combined-accessions.txt"

TEMP="${ACCESSION_DIR}/all-accessions.tmp"

# Get header from first BioProject

head -n 1 \
"${SEARCH_DIR}/${PROJECTS[0]}/bactopia-accessions.txt" \
> "$COMBINED"

# Combine accession records

for PROJECT in "${PROJECTS[@]}"; do

    tail -n +2 \
    "${SEARCH_DIR}/${PROJECT}/bactopia-accessions.txt" \
    >> "$TEMP"

done

# Remove duplicate Experiment accessions

awk -F'\t' '!seen[$1]++' "$TEMP" >> "$COMBINED"

rm -f "$TEMP"

# ============================================================
# 7. Count unique samples
# ============================================================

TOTAL=$(awk 'END {print NR-1}' "$COMBINED")

echo "=========================================="
echo "Total unique Experiment accessions: $TOTAL"
echo "=========================================="

if [[ "$TOTAL" -eq 0 ]]; then

    echo "ERROR: No samples retrieved."

    exit 1

fi

# ============================================================
# 8. Run Bactopia
# ============================================================

echo "Starting Bactopia analysis..."

cd "$BASE_DIR"

if bactopia \
    --accessions "$COMBINED" \
    --coverage 100 \
    --max_cpus 8 \
    --outdir "$RESULT_DIR" \
    -w "$WORK_DIR"

then

    echo "=========================================="
    echo "Bactopia completed successfully!"
    echo "=========================================="

else

    echo "=========================================="
    echo "Bactopia failed!"
    echo "=========================================="

    exit 1

fi

# ============================================================
# 9. Final summary
# ============================================================

echo "Total unique samples: $TOTAL"

echo "Results directory:"
echo "$RESULT_DIR"

echo "Workflow completed."