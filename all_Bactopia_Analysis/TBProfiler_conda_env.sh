#!/bin/bash
#SBATCH --job-name=TBprofiler_Conda       # Job name
#SBATCH --partition=batch                  # Queue/partition
#SBATCH --ntasks=1                         # Number of tasks
#SBATCH --cpus-per-task=8                  # Number of CPU cores
#SBATCH --mem=40gb                         # Memory allocation
#SBATCH --time=07-00:00:00                 # Max runtime (7 days)
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out   # STDOUT log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err    # STDERR log
#SBATCH --mail-type=END,FAIL               # Email notifications
#SBATCH --mail-user=ma95362@uga.edu       # Email recipient

set -euo pipefail  # safer bash settings

# Initialize Micromamba (adjust path if needed)
source /apps/eb/Micromamba/2.3.0/etc/profile.d/conda.sh
micromamba activate tbprofiler

# Set working directories
FASTQ_DIR="/scratch/ma95362/musse_MGA/merged/first_run_merged/all_reads"
OUTDIR="/scratch/ma95362/my_tbprofiler_results"
mkdir -p "$OUTDIR"

# Log file for processed samples
LOGFILE="$OUTDIR/processed_samples.log"
echo "TB-Profiler run started at $(date)" > "$LOGFILE"

# Loop through all R1 FASTQ files
for read1 in "$FASTQ_DIR"/*._R1.fastq.gz; do
    # Get base sample name
    base=$(basename "$read1" ._R1.fastq.gz)
    read2="$FASTQ_DIR/${base}._R2.fastq.gz"
    
    # Check that R2 exists
    if [[ ! -f "$read2" ]]; then
        echo "WARNING: $read2 not found, skipping $base" >> "$LOGFILE"
        continue
    fi
    
    SAMPLE_OUT="$OUTDIR/$base"
    mkdir -p "$SAMPLE_OUT"
    
    echo "Processing sample: $base" | tee -a "$LOGFILE"
    
    tb-profiler profile \
        -1 "$read1" \
        -2 "$read2" \
        -p "$base" \
        --db /home/ma95362/.conda/envs/mtbvartools/share/tbprofiler/tbdb \
        --spoligotype \
        --dir "$SAMPLE_OUT" \
        --threads 8 \
        --txt
done

echo "TB-Profiler run finished at $(date)" >> "$LOGFILE"

# Deactivate the environment
micromamba deactivate

