#!/bin/bash
#SBATCH --job-name=Abebe_Tbprofiler_Lineages_Spoligo         # Job name
#SBATCH --partition=batch                                    # Partition (queue) name
#SBATCH --ntasks=1                                           # Run on a single CPU
#SBATCH --cpus-per-task=8                                    # Number of cores per task
#SBATCH --mem=40gb                                           # Job memory request
#SBATCH --time=05-00:00:00                                   # Time limit d-hh:mm:ss
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out         # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err          # Standard error log
#SBATCH --mail-type=END,FAIL                                 # Mail notifications
#SBATCH --mail-user=ma95362@uga.edu                          # Your email address

# Load TBProfiler
module load TBProfiler/6.6.2

#move to working directory
cd $OUTDIR

# Define input and output directories
BASE_DIR="/scratch/ma95362/ETH_bovis_Sequence/Abebe_all_samples/snippy_results"
OUT_DIR="/scratch/ma95362/ETH_bovis_Sequence/Abebe_all_samples/tbprofiler_results"
DB_DIR="/scratch/ma95362/ETH_bovis_Sequence/tbprofiler_db"  # 🛠️ Local db path (writable)

# Create necessary dirs
mkdir -p "$OUT_DIR"
mkdir -p "$DB_DIR"

# Initialize DB (only once)
if [ ! -f "$DB_DIR/tbdb.js" ]; then
    echo "[INFO] Creating local TBProfiler DB at $DB_DIR"
    tb-profiler update_tbdb --db-dir "$DB_DIR"
fi

# Run TBProfiler for each BAM
for SAMPLE_DIR in "$BASE_DIR"/*; do
    SAMPLE=$(basename "$SAMPLE_DIR")

    if [ "$SAMPLE" == "bactopia-runs" ] || [ ! -d "$SAMPLE_DIR" ]; then
        continue
    fi

    BAM="$SAMPLE_DIR/tools/snippy/genomic/${SAMPLE}.bam"
    if [ -f "$BAM" ]; then
        echo "Running TBProfiler for $SAMPLE"
        tb-profiler profile \
            --bam "$BAM" \
            --prefix "$OUT_DIR/${SAMPLE}" \
            --db-dir "$DB_DIR"
    else
        echo "⚠️ BAM file not found for $SAMPLE — skipping"
    fi
done