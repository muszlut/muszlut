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


# Define input and output directories
BASE_DIR="/scratch/ma95362/ETH_bovis_Sequence/Abebe_all_samples/snippy_results"
OUT_DIR="/scratch/ma95362/ETH_bovis_Sequence/Abebe_all_samples/tbprofiler_results"
mkdir -p "$OUT_DIR"

# Load TBProfiler
module load TBProfiler/6.6.2

#move to working directory
cd $OUTDIR

# Loop through each sample
for SAMPLE_DIR in "$BASE_DIR"/*; do
    SAMPLE=$(basename "$SAMPLE_DIR")

    # Skip non-directories (e.g., bactopia-runs)
    if [ ! -d "$SAMPLE_DIR" ] || [ "$SAMPLE" == "bactopia-runs" ]; then
        continue
    fi

    BAM="$SAMPLE_DIR/tools/snippy/genomic/${SAMPLE}.bam"
    if [ -f "$BAM" ]; then
        echo "Running TBProfiler for $SAMPLE"
        tb-profiler profile \
            --bam "$BAM" \
            --prefix "$OUT_DIR/${SAMPLE}" \
            --force
    else
        echo "⚠️ BAM file not found for $SAMPLE — skipping"
    fi
done

# Optional: Collate all sample results into a summary
tb-profiler collate \
    --dir "$OUTDIR" \
    --prefix "$OUTDIR/tbprofiler_summary"
