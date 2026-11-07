#!/bin/bash
#SBATCH --job-name=Mtb_PRJNA1174701
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/ma95362/Bio_project_publication/Mtb_Multiple_BioProjects/logs/PRJNA1174701.%j.out
#SBATCH --error=/scratch/ma95362/Bio_project_publication/Mtb_Multiple_BioProjects/logs/PRJNA1174701.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# PRJNA1174701: Bactopia's downloader fails for every sample ("not found"),
# but ENA serves the FASTQs. So:
#   1. download R1/R2 directly from ENA (skip samples already in results_singularity)
#   2. name them by Experiment accession (SRX...), like the other samples
#   3. process them as local FASTQs into results_singularity (same settings as main run)

set -euo pipefail

module load Mamba/23.11.0-0
source activate bactopia4
export PYTHONNOUSERSITE=1
module load Apptainer 2>/dev/null || module load Singularity 2>/dev/null || true
export NXF_SINGULARITY_CACHEDIR=/scratch/ma95362/singularity_cache

BASE_DIR="/scratch/ma95362/Bio_project_publication/Mtb_Multiple_BioProjects"
PROJECT="PRJNA1174701"
RUN_DIR="${BASE_DIR}/runs/${PROJECT}"        # own launch folder / Nextflow cache
FASTQ_DIR="${RUN_DIR}/fastqs"
WORK_DIR="${BASE_DIR}/work_${PROJECT}"
RESULT_DIR="${BASE_DIR}/results_singularity"
SAMPLES="${RUN_DIR}/${PROJECT}-samples.txt"

mkdir -p "$RUN_DIR" "$FASTQ_DIR" "$WORK_DIR" "${BASE_DIR}/logs"
cd "$RUN_DIR"

# ------------------------------------------------------------
# 1. File list from ENA: run, experiment, FASTQ URLs, md5s
# ------------------------------------------------------------
REPORT="${RUN_DIR}/${PROJECT}-ena-filereport.tsv"
curl -s --retry 5 \
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${PROJECT}&result=read_run&fields=run_accession,experiment_accession,library_layout,fastq_ftp,fastq_md5&format=tsv" \
  > "$REPORT"
echo "Runs listed by ENA: $(( $(wc -l < "$REPORT") - 1 ))"

# ------------------------------------------------------------
# 2. Download R1/R2 per run (skip unpaired 3rd file), verify md5,
#    then join per Experiment as SRXxxx_R1/_R2.fastq.gz
# ------------------------------------------------------------
tail -n +2 "$REPORT" | while IFS=$'\t' read -r RUN EXP LAYOUT FTP MD5; do

    # already processed in the main run?
    [[ -d "${RESULT_DIR}/${EXP}" ]] && continue
    # already downloaded (from an earlier attempt)?
    [[ -s "${FASTQ_DIR}/${EXP}_R1.fastq.gz" && -f "${FASTQ_DIR}/.${RUN}.done" ]] && continue

    IFS=';' read -ra URLS <<< "$FTP"
    IFS=';' read -ra SUMS <<< "$MD5"
    for i in "${!URLS[@]}"; do
        U="${URLS[$i]}"
        case "$U" in
            *_1.fastq.gz) R=R1 ;;
            *_2.fastq.gz) R=R2 ;;
            *) continue ;;           # skip unpaired SRRxxx.fastq.gz
        esac
        TMP="${FASTQ_DIR}/${RUN}_${R}.fastq.gz"
        for TRY in 1 2 3 4 5; do
            wget -q -O "$TMP" "https://${U}" && \
              echo "${SUMS[$i]}  ${TMP}" | md5sum -c --status && break
            echo "Retry ${TRY} for ${RUN} ${R}"; sleep 30
        done
        echo "${SUMS[$i]}  ${TMP}" | md5sum -c --status || { echo "FAILED: ${RUN} ${R}"; rm -f "$TMP"; continue 2; }
        cat "$TMP" >> "${FASTQ_DIR}/${EXP}_${R}.fastq.gz"   # gz files can be concatenated
        rm -f "$TMP"
    done
    touch "${FASTQ_DIR}/.${RUN}.done"
    echo "Downloaded: ${EXP} (${RUN})"
done

N_DL=$(ls "${FASTQ_DIR}"/*_R1.fastq.gz 2>/dev/null | wc -l)
echo "=========================================="
echo "Experiments with FASTQs ready: ${N_DL}"
echo "=========================================="
[[ "$N_DL" -gt 0 ]] || { echo "Nothing to process."; exit 0; }

# ------------------------------------------------------------
# 3. Sample sheet from the downloaded FASTQs
# ------------------------------------------------------------
bactopia prepare \
    --path "$FASTQ_DIR" \
    --fastq-ext .fastq.gz \
    --species "Mycobacterium tuberculosis" \
    --genome-size 4411532 \
    > "$SAMPLES"

echo "Samples in sheet: $(( $(wc -l < "$SAMPLES") - 1 ))"
tail -n +2 "$SAMPLES" | cut -f2 | sort | uniq -c

# ------------------------------------------------------------
# 4. Run Bactopia (same settings as the main run)
# ------------------------------------------------------------
CUSTOM_CONFIG="${RUN_DIR}/custom.config"
cat > "$CUSTOM_CONFIG" <<'EOF'
process {
    withName: '.*' {
        time          = '12.h'
        errorStrategy = { task.attempt <= 1 ? 'retry' : 'ignore' }
        maxRetries    = 1
    }
}
EOF

bactopia \
    --samples "$SAMPLES" \
    --coverage 100 \
    --max_cpus 8 \
    --max_memory 100.GB \
    --max_time 12.h \
    --outdir "$RESULT_DIR" \
    -w "$WORK_DIR" \
    -c "$CUSTOM_CONFIG" \
    -profile singularity \
    -resume

echo "Done. PRJNA1174701 results are in ${RESULT_DIR}"