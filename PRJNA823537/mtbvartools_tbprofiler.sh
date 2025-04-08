#!/bin/bash
#SBATCH --job-name=mtbvarTools_tbprofiler_job         # Job name
#SBATCH --partition=batch                             # Partition (queue) name
#SBATCH --ntasks=1                                    # Run on a single CPU
#SBATCH --cpus-per-task=8                             # Number of cores per task
#SBATCH --mem=40gb                                    # Job memory request
#SBATCH --time=07-00:00:00                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out  # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err   # Standard error log

#SBATCH --mail-type=END,FAIL                          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                   # Where to send mail

# Load modules or set up environment
module load Python/3.9.6-GCCcore-11.2.0               # Adjust Python version based on your setup
source ~/mtbvartools_env/bin/activate                 # Activate MTBvarTools environment

# Define input/output directories
BAM_DIR="/scratch/ma95362/bam_files"                  # Directory containing BAM files
OUT_DIR="/scratch/ma95362/mtbvartools_results"        # Directory for results

# Create output directory if it doesn't exist
if [ ! -d "$OUT_DIR" ]; then
    mkdir -p "$OUT_DIR"
fi

# Loop through BAM files and run MTBvarTools profile
for BAM_FILE in "$BAM_DIR"/*.bam; do
    SAMPLE_NAME=$(basename "$BAM_FILE" .bam)
    echo "Processing sample: $SAMPLE_NAME"
    
    # Run MTBvarTools profile for lineage prediction and drug resistance
    mtbvarTools profile --bam "$BAM_FILE" --out "$OUT_DIR/$SAMPLE_NAME"
    if [ $? -ne 0 ]; then
        echo "Error processing $SAMPLE_NAME. Skipping..."
        continue
    fi

    echo "Finished processing sample: $SAMPLE_NAME"
done

echo "All samples processed successfully!"
