#!/bin/bash
#SBATCH --job-name=fixing_FASTQs                   # Job name
#SBATCH --partition=batch                           # Partition (queue)
#SBATCH --ntasks=1                                  # Single task
#SBATCH --cpus-per-task=8                           # CPUs per task
#SBATCH --mem=40gb                                  # Memory
#SBATCH --time=05-00:00:00                          # Time limit (HH:MM:SS)
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out    # STDOUT log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err     # STDERR log
#SBATCH --mail-type=END,FAIL                         # Mail events
#SBATCH --mail-user=ma95362@uga.edu                  # Your email

#-----------------------------------
# Modules
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
# Log file
#-----------------------------------
LOGFILE="$FIXED_DIR/fastq_fix_report.txt"
echo "FASTQ Fix Report - $(date)" > "$LOGFILE"
echo "----------------------------------------" >> "$LOGFILE"

#-----------------------------------
# Function to fix a FASTQ and log discarded reads
#-----------------------------------
fix_fastq() {
    local infile=$1
    local outfile=$2

    if [[ ! -s "$infile" ]]; then
        echo "$infile is empty, skipping" | tee -a "$LOGFILE"
        return
    fi

    echo "Processing $infile ..."
    reformat.sh in="$infile" out="$outfile" tossbrokenreads 2> tmp.log

    # Extract number of discarded reads from stderr
    if grep -q "Discarded" tmp.log; then
        discarded=$(grep "Discarded" tmp.log | awk '{print $2}')
    else
        discarded=0
    fi

    echo "$infile -> $outfile : $discarded reads discarded" | tee -a "$LOGFILE"
    rm -f tmp.log
}

#-----------------------------------
# Process paired-end FASTQs
#-----------------------------------
for fq1 in *_R1.fastq.gz; do
    fq2="${fq1/_R1/_R2}"  # Corresponding R2
    if [[ -f "$fq2" ]]; then
        fix_fastq "$fq1" "$FIXED_DIR/$fq1"
        fix_fastq "$fq2" "$FIXED_DIR/$fq2"
    else
        # Single-end read
        fix_fastq "$fq1" "$FIXED_DIR/$fq1"
    fi
done

#-----------------------------------
# Process remaining single-end FASTQs not caught above
#-----------------------------------
for fq in *_R2.fastq.gz; do
    if [[ ! -f "${fq/_R2/_R1}" ]]; then
        fix_fastq "$fq" "$FIXED_DIR/$fq"
    fi
done

#-----------------------------------
echo "All FASTQs processed. Fixed files are in $FIXED_DIR"
echo "Detailed report saved in $LOGFILE"
