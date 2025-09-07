#!/bin/bash
#SBATCH --job-name=fixing_FASTQs
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40gb
#SBATCH --time=05-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

#-----------------------------------
# Load BBMap
#-----------------------------------
module load BBMap/39.19

#-----------------------------------
# Directories
#-----------------------------------
WORKDIR="/scratch/ma95362/ETH_bovis_Sequence/reads/64_SH_Sequence_data/raw"
FIXED_DIR="$WORKDIR/fixed_fastqs"
mkdir -p "$FIXED_DIR"
cd "$WORKDIR" || { echo "Cannot cd to $WORKDIR"; exit 1; }

#-----------------------------------
# Log files
#-----------------------------------
SUMMARY="$FIXED_DIR/fastq_summary_report.txt"
DETAIL="$FIXED_DIR/fastq_fix_report.txt"
REPAIR="$FIXED_DIR/fastq_repair_report.txt"

echo "FASTQ Summary Report - $(date)" > "$SUMMARY"
echo "----------------------------------------" >> "$SUMMARY"
echo "FASTQ Fix Report - $(date)" > "$DETAIL"
echo "----------------------------------------" >> "$DETAIL"
echo "FASTQ Repair Report - $(date)" > "$REPAIR"
echo "----------------------------------------" >> "$REPAIR"

#-----------------------------------
# Function to process paired FASTQs
#-----------------------------------
process_pair() {
    local fq1=$1
    local fq2=$2
    local base=$(basename "$fq1" _R1.fastq.gz)

    local out1="$FIXED_DIR/${base}_R1.fastq.gz"
    local out2="$FIXED_DIR/${base}_R2.fastq.gz"
    local singles="$FIXED_DIR/${base}_singletons.fastq.gz"

    echo "Processing $fq1 and $fq2 ..." | tee -a "$DETAIL"

    # Run reformat.sh and capture log
    tmp_log=$(mktemp)
    reformat.sh in="$fq1" in2="$fq2" out1="$out1" out2="$out2" tossbrokenreads overwrite 2> "$tmp_log"

    # Extract discarded reads
    discarded=$(grep -oP 'Discarded reads: \K\d+' "$tmp_log" || echo 0)
    echo "$base : $discarded reads discarded" | tee -a "$DETAIL"

    # If any reads discarded, run repair.sh
    if [[ "$discarded" -gt 0 ]]; then
        echo "Pair mismatch detected for $base → running repair.sh" | tee -a "$REPAIR"
        repair.sh in1="$fq1" in2="$fq2" out1="$out1" out2="$out2" outs="$singles" overwrite 2>> "$REPAIR"
    fi

    rm -f "$tmp_log"
}

#-----------------------------------
# Main loop
#-----------------------------------
total_pairs=0
for fq1 in *_R1.fastq.gz; do
    fq2="${fq1/_R1/_R2}"
    if [[ -f "$fq2" ]]; then
        process_pair "$fq1" "$fq2"
        ((total_pairs++))
    fi
done

echo "Total pairs processed: $total_pairs" >> "$SUMMARY"
echo "All FASTQs processed. Reports saved in $FIXED_DIR"
