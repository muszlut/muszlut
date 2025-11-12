#!/bin/bash
#SBATCH --job-name=Split-and_copy_reads
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120gb
#SBATCH --time=03-00:00:00
#SBATCH --output=/scratch/ma95362/EPTB_Hilina/Newe/logs/log.%j.out
#SBATCH --error=/scratch/ma95362/EPTB_Hilina/Newe/logs/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Optional: Load BBTools if needed
# module load BBMap
module load BBMap/39.19-GCC-13.3.0
# Set base directory
BASE_DIR="/scratch/ma95362/EPTB_Hilina/Bactopia_Run"
DEST_DIR="/scratch/ma95362/EPTB_Hilina/Bactopia_Run/split_reads"
cd "$BASE_DIR"
# Create destination directory if it doesn't exist
mkdir -p "$DEST_DIR"
# Loop through all SRX sample folders
for sample in "$BASE_DIR"/SRX*/; do
  sample_id=$(basename "$sample")
  fq="$sample/main/qc/${sample_id}.fastq.gz"

  if [[ -f "$fq" ]]; then
    echo "🔧 Splitting $sample_id..."
    reformat.sh in="$fq" out1="${sample_id}_R1.fastq.gz" out2="${sample_id}_R2.fastq.gz"

    echo "📦 Copying to $DEST_DIR..."
    mv "${sample_id}_R1.fastq.gz" "$DEST_DIR/"
    mv "${sample_id}_R2.fastq.gz" "$DEST_DIR/"
  else
    echo "⚠️ FASTQ not found for $sample_id"
  fi
done

echo "✅ All samples processed and copied to $DEST_DIR"