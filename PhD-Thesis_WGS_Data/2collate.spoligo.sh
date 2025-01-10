#!/bin/bash
#SBATCH --job-name=Spoligo_collate       # Job name
#SBATCH --partition=batch                # Partition (queue) name
#SBATCH --ntasks=1                       # Run on a single CPU
#SBATCH --cpus-per-task=8                # Number of cores per task
#SBATCH --mem=40gb                       # Job memory request
#SBATCH --time=06:00:00                  # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out     # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err      # Standard error log

#SBATCH --mail-type=END,FAIL             # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu      # Where to send mail

# Load jq module (if available on your system)
module load jq

# Define the directory containing the spoligotype JSON files
JSON_DIR="/scratch/ma95362/musse_MGA/fastqs/MGA_paired_end_samples/results"
# Define the output CSV file
OUTPUT_CSV="/scratch/ma95362/musse_MGA/fastqs/MGA_paired_end_samples/results/collated_spoligotype.csv"

# Add headers to the CSV file
echo "SampleName,Binary,Octal,Family,SIT,Countries,SpacerName,SpacerSeq,SpacerCount" > "$OUTPUT_CSV"

# Loop through each JSON file in the directory
for json_file in "$JSON_DIR"/*.spoligotype.json; do
    # Extract the sample name from the file name
    SAMPLE_NAME=$(basename "$json_file" .spoligotype.json)
    
    # Use jq to parse the JSON and format as CSV
    jq -r --arg sample "$SAMPLE_NAME" '
    [
        $sample,
        (.binary | @sh),
        .octal,
        .family,
        .SIT,
        .countries
    ] as $common |
    .spacers[] | [$common[], .name, .seq, .count] | @csv' "$json_file" >> "$OUTPUT_CSV"
done

# Print a message indicating completion
echo "Collation complete. The collated CSV file is located at $OUTPUT_CSV"
