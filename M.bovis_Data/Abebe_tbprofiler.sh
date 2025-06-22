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


# Input and output directories
SNIPPY_DIR="/scratch/ma95362/ETH_bovis_Sequence/Abebe_all_samples/snippy_results"
OUTDIR="/scratch/ma95362/ETH_bovis_Sequence/Abebe_all_samples/tbprofiler_results"

mkdir -p "$OUTDIR"

# Load TBProfiler
module load TBProfiler/6.6.2

#move to working directory
cd $OUTDIR

# Loop through each Snippy result directory
for SAMPLE in "$SNIPPY_DIR"/*; do
    if [ -d "$SAMPLE" ]; then
        SAMPLE_ID=$(basename "$SAMPLE")
        BAM="$SAMPLE/snps.bam"

        if [ -f "$BAM" ]; then
            echo "Processing $SAMPLE_ID..."
            tb-profiler profile \
                --bam "$BAM" \
                --prefix "$OUTDIR/$SAMPLE_ID" \
                --threads 8 \
                --force
        else
            echo "⚠️ BAM not found for $SAMPLE_ID — skipping"
        fi
    fi
done

# Optional: Collate all sample results into a summary
tb-profiler collate \
    --dir "$OUTDIR" \
    --prefix "$OUTDIR/tbprofiler_summary"
