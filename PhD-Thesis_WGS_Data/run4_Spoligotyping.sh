#!/bin/bash
#SBATCH --job-name=Musse_WGS_spoligo_collate           # Job name
#SBATCH --partition=batch                      # Partition (queue) name
#SBATCH --ntasks=1                             # Run on a single CPU
#SBATCH --cpus-per-task=6                      # Number of cores per task
#SBATCH --mem=10gb                             # Job memory request
#SBATCH --time=12:00:00                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out   # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu            # Where to send mail

# Set output directory variable
OUTDIR="/scratch/ma95362/musse_MGA/fastqs/MGA_paired_end_samples"
# Define path to the sequence reads
FASTQ_DIR="/work/fdqlab/Ethiopia_wgs_mtb_2024/first_run"

# Create output directory if it doesn't exist
if [ ! -d "$OUTDIR" ]; then
    mkdir -p "$OUTDIR"
fi

# Load necessary modules
module load Miniconda3/23.5.2-0
module load jq

# Activate conda environment (should be version 6.3.0)
source activate tb-profiler-env

# Move to working directory
cd "$OUTDIR"

# Initialize output files
COLLATED_CSV="$OUTDIR/collated_spoligotype.csv"
COLLATED_TXT="$OUTDIR/collated_spoligotype.txt"

# Add headers to the collated CSV file
echo "SampleName,Key,Value" > "$COLLATED_CSV"

# Loop through all pairs of read files and run tb-profiler spoligotype analysis
for READ1 in "$FASTQ_DIR"/*_R1.fastq.gz; do
    # Generate the corresponding READ2 file name
    READ2="${READ1/_R1.fastq.gz/_R2.fastq.gz}"
    
    # Extract sample name from file name
    SAMPLE_NAME=$(basename "$READ1" _R1.fastq.gz)
    
    # Run tb-profiler spoligotype analysis with paired-end reads
    tb-profiler spoligotype --read1 "$READ1" --read2 "$READ2" --prefix "$SAMPLE_NAME"
    
    # Convert JSON to CSV and append to collated CSV file
    jq -r --arg sample "$SAMPLE_NAME" 'to_entries | map([$sample, .key, (.value | tostring)] | join(",")) | join("\n")' "$OUTDIR/$SAMPLE_NAME.spoligotype.json" >> "$COLLATED_CSV"
    
    # Convert JSON to TXT and append to collated TXT file
    echo "Sample: $SAMPLE_NAME" >> "$COLLATED_TXT"
    jq -r '
        "Binary: " + .binary + "\n" +
        "Octal: " + .octal + "\n" +
        "Family: " + .family + "\n" +
        "SIT: " + .SIT + "\n" +
        "Countries: " + .countries + "\n\n" +
        "Spacers:\n" +
        ( .spacers[] | "Name: " + .name + "\nSeq: " + .seq + "\nCount: " + (.count | tostring) + "\n" )
    ' "$OUTDIR/$SAMPLE_NAME.spoligotype.json" >> "$COLLATED_TXT"
    echo "" >> "$COLLATED_TXT"
done

# Print a message indicating completion
echo "Spoligotyping analysis and JSON-to-TXT/CSV conversion complete. Collated results are in $OUTDIR"
