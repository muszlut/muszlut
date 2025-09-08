#!/bin/bash
#SBATCH --job-name=kraken2_summary
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Load modules
module load Kraken2/2.1.3-gompi-2023a
module load Python/3.11.5

# Path to your KrakenTools installation
KRAKENTOOLS=/scratch/ma95362/KrakenTools   # <-- update if needed
REPORTS_DIR=/scratch/ma95362/ETH_bovis_Sequence/reads/64_SH_Sequence_data/raw/fixed_fastqs/fixed_fastqs
OUTFILE=/scratch/ma95362/ETH_bovis_Sequence/reads/64_SH_Sequence_data/raw/fixed_fastqs/fixed_fastqs/kraken2_combined.tsv

# Collect all report files
# -------------------------------
REPORT_FILES=(${REPORTS_DIR}/*_kraken2_report.txt)

# Generate sample names (strip path and suffix)
SAMPLE_NAMES=()
for f in "${REPORT_FILES[@]}"; do
    name=$(basename "$f" _kraken2_report.txt)
    SAMPLE_NAMES+=("$name")
done

# -------------------------------
# Run KrakenTools to combine reports
# -------------------------------
python $KRAKENTOOLS/combine_kreports.py \
    -r "${REPORT_FILES[@]}" \
    -o "$OUTFILE" \
    --sample-names "${SAMPLE_NAMES[@]}"

echo "✅ Combined Kraken2 report written to: $OUTFILE"