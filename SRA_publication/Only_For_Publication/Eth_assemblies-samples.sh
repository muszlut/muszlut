#!/bin/bash
#SBATCH --job-name=Mtb_assembly_PRJNA975069
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/ma95362/Bio_project_publication/Mtb_Multiple_BioProjects/logs/PRJNA975069.%j.out
#SBATCH --error=/scratch/ma95362/Bio_project_publication/Mtb_Multiple_BioProjects/logs/PRJNA975069.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# PRJNA975069 has assemblies only (no SRA reads).
# This script: downloads the assemblies from NCBI -> builds a Bactopia
# sample sheet (runtype = assembly) -> runs Bactopia, which simulates
# reads from each assembly and processes them like the other samples.

set -euo pipefail

# ============================================================
# 1. Environment
# ============================================================
module load Mamba/23.11.0-0
source activate bactopia4
export PYTHONNOUSERSITE=1

module load Apptainer 2>/dev/null || module load Singularity 2>/dev/null || true
export NXF_SINGULARITY_CACHEDIR=/scratch/ma95362/singularity_cache
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

# NCBI datasets CLI is needed for the download.
# If this check fails, install it once (on a login node):
#   source activate bactopia4
#   mamba install -c conda-forge ncbi-datasets-cli
if ! command -v datasets >/dev/null 2>&1; then
    module load ncbi-datasets-cli 2>/dev/null || true
fi
command -v datasets >/dev/null 2>&1 || { echo "ERROR: 'datasets' not found. See install note in script."; exit 1; }

# ============================================================
# 2. Directories and settings
# ============================================================
BASE_DIR="/scratch/ma95362/Bio_project_publication/Mtb_Multiple_BioProjects"
PROJECT="PRJNA975069"

ASM_DIR="${BASE_DIR}/assemblies/${PROJECT}"
FNA_DIR="${ASM_DIR}/fna"
SAMPLES="${ASM_DIR}/${PROJECT}-samples.txt"

RESULT_DIR="${BASE_DIR}/results"        # same results folder as the other BioProjects
WORK_DIR="${BASE_DIR}/work_${PROJECT}"

GENOME_SIZE=4411532
SPECIES="Mycobacterium tuberculosis"

mkdir -p "$ASM_DIR" "$FNA_DIR" "$RESULT_DIR" "$WORK_DIR" "${BASE_DIR}/logs"

# ============================================================
# 3. Download all assemblies linked to the BioProject
# ============================================================
cd "$ASM_DIR"

if [[ ! -s "${PROJECT}.zip" ]]; then
    echo "Downloading assemblies for ${PROJECT} ..."
    datasets download genome accession "$PROJECT" \
        --assembly-source GenBank \
        --include genome \
        --filename "${PROJECT}.zip"
fi

rm -rf ncbi_dataset
unzip -q -o "${PROJECT}.zip"

# ============================================================
# 4. One gzipped FASTA per sample, named by assembly accession
#    e.g. GCA_030415455.1.fna.gz
# ============================================================
echo "Collecting FASTA files ..."

for DIR in ncbi_dataset/data/GCA_*; do
    ACC=$(basename "$DIR")
    FNA=$(ls "$DIR"/*.fna 2>/dev/null | head -n 1 || true)
    if [[ -n "$FNA" ]]; then
        gzip -c "$FNA" > "${FNA_DIR}/${ACC}.fna.gz"
    else
        echo "WARNING: no FASTA for ${ACC}"
    fi
done

N_ASM=$(ls "${FNA_DIR}"/*.fna.gz | wc -l)
echo "=========================================="
echo "Assemblies collected: ${N_ASM} (BioProject lists 146 BioSamples)"
echo "=========================================="

# Lookup table: assembly accession -> BioSample -> strain name (for your metadata)
if command -v dataformat >/dev/null 2>&1; then
    dataformat tsv genome \
        --inputfile ncbi_dataset/data/assembly_data_report.jsonl \
        --fields accession,assminfo-biosample-accession,organism-infraspecific-strain \
        > "${ASM_DIR}/${PROJECT}-accession-to-strain.tsv" || true
fi

# ============================================================
# 5. Build the Bactopia sample sheet
# ============================================================
echo "Preparing sample sheet ..."

bactopia prepare \
    --path "$FNA_DIR" \
    --assembly-ext .fna.gz \
    --species "$SPECIES" \
    --genome-size "$GENOME_SIZE" \
    > "$SAMPLES"

TOTAL=$(( $(wc -l < "$SAMPLES") - 1 ))
echo "Samples in sheet: $TOTAL"
echo "Run types:"
tail -n +2 "$SAMPLES" | cut -f2 | sort | uniq -c

if [[ "$TOTAL" -le 0 ]]; then
    echo "ERROR: empty sample sheet."
    exit 1
fi

# ============================================================
# 6. Run Bactopia on the assemblies
# ============================================================
cd "$BASE_DIR"

bactopia \
    --samples "$SAMPLES" \
    --max_cpus 32 \
    --max_memory 120.GB \
    --outdir "$RESULT_DIR" \
    -w "$WORK_DIR" \
    -profile singularity \
    -resume

echo "=========================================="
echo "Done: ${TOTAL} PRJNA975069 assemblies processed"
echo "Results: ${RESULT_DIR}"
echo "=========================================="