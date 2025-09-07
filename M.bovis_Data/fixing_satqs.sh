#!/bin/bash
#SBATCH --job-name=check_repair_FASTQs       # Job name
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
DRYRUN_LOG="$FIXED_DIR/fastq_dryrun_report.txt"
REPAIR_LOG="$FIXED_DIR/fastq_repair_report.txt"
SUMMARY_LOG="$FIXED_DIR/fastq_summary_report.txt"

echo "FASTQ Dry-Run Report - $(date)" > "$DRYRUN_LOG"
echo "----------------------------------------" >> "$DRYRUN_LOG"
echo "FASTQ Repair Report - $(date)" > "$REPAIR_LOG"
echo "----------------------------------------" >> "$REPAIR_LOG"
echo "FASTQ Summary Report - $(date)" > "$SUMMARY_LOG"
echo "----------------------------------------" >> "$SUMMARY_LOG"

#-----------------------------------
# Function to check broken reads
#-----------------------------------
check_broken_fastq() {
    local fq1=$1
    local fq2=$2
    local broken=0

    if [[ -z "$fq2" ]]; then
        reformat.sh in="$fq1" out=/dev/null tossbrokenreads 2> tmp.log
    else
        reformat.sh in="$fq1" in2="$fq2" out=/dev/null out2=/dev/null repair 2> tmp.log
    fi

    if grep -q "Discarded" tmp.log; then
        broken=$(grep "Discarded" tmp.log | awk '{print $2}')
    fi
    rm -f tmp.log
    echo $broken
}

#-----------------------------------
# Keep track of repaired files
#-----------------------------------
repaired_files=0
total_discarded=0

#-----------------------------------
# Step 1 & 2: Dry-run + repair only if broken reads found
#-----------------------------------
for fq1 in *_R1.fastq.gz; do
    fq2="${fq1/_R1/_R2}"
    if [[ -f "$fq2" ]]; then
        broken=$(check_broken_fastq "$fq1" "$fq2")
        echo "$fq1 & $fq2 -> Broken reads: $broken" | tee -a "$DRYRUN_LOG"
        if [[ $broken -gt 0 ]]; then
            echo "Repairing $fq1 & $fq2 ..."
            out1="$FIXED_DIR/${fq1%.fastq.gz}_fixed.fastq.gz"
            out2="$FIXED_DIR/${fq2%.fastq.gz}_fixed.fastq.gz"
            reformat.sh in="$fq1" in2="$fq2" out="$out1" out2="$out2" repair overwrite 2> tmp.log
            repaired=$(grep "Discarded" tmp.log | awk '{print $2}')
            echo "$fq1 & $fq2 -> Repaired, discarded: $repaired" | tee -a "$REPAIR_LOG"
            ((repaired_files++))
            ((total_discarded+=repaired))
            rm -f tmp.log
        fi
    else
        broken=$(check_broken_fastq "$fq1")
        echo "$fq1 -> Broken reads: $broken" | tee -a "$DRYRUN_LOG"
        if [[ $broken -gt 0 ]]; then
            echo "Repairing $fq1 ..."
            out1="$FIXED_DIR/${fq1%.fastq.gz}_fixed.fastq.gz"
            reformat.sh in="$fq1" out="$out1" tossbrokenreads overwrite 2> tmp.log
            repaired=$(grep "Discarded" tmp.log | awk '{print $2}')
            echo "$fq1 -> Repaired, discarded: $repaired" | tee -a "$REPAIR_LOG"
            ((repaired_files++))
            ((total_discarded+=repaired))
            rm -f tmp.log
        fi
    fi
done

#-----------------------------------
# Handle remaining single-end R2 files without R1
#-----------------------------------
for fq in *_R2.fastq.gz; do
    if [[ ! -f "${fq/_R2/_R1}" ]]; then
        broken=$(check_broken_fastq "$fq")
        echo "$fq -> Broken reads: $broken" | tee -a "$DRYRUN_LOG"
        if [[ $broken -gt 0 ]]; then
            echo "Repairing $fq ..."
            out="$FIXED_DIR/${fq%.fastq.gz}_fixed.fastq.gz"
            reformat.sh in="$fq" out="$out" tossbrokenreads overwrite 2> tmp.log
            repaired=$(grep "Discarded" tmp.log | awk '{print $2}')
            echo "$fq -> Repaired, discarded: $repaired" | tee -a "$REPAIR_LOG"
            ((repaired_files++))
            ((total_discarded+=repaired))
            rm -f tmp.log
        fi
    fi
done

#-----------------------------------
# Summary
#-----------------------------------
echo "Total repaired files: $repaired_files" | tee -a "$SUMMARY_LOG"
echo "Total discarded reads: $total_discarded" | tee -a "$SUMMARY_LOG"
echo "All done. Repaired FASTQs (if any) are in $FIXED_DIR"
echo "Dry-run log: $DRYRUN_LOG"
echo "Repair log: $REPAIR_LOG"
echo "Summary log: $SUMMARY_LOG"
