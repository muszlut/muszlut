#!/bin/bash
#SBATCH --job-name=tb_profile_parsing                          # Job name
#SBATCH --partition=batch                                      # Partition (queue) name
#SBATCH --ntasks=1                                             # Run on a single CPU
#SBATCH --cpus-per-task=8                                      # Number of cores per task
#SBATCH --mem=40gb                                             # Job memory request
#SBATCH --time=07-00:00:00                                     # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err            # Standard error log

#SBATCH --mail-type=END,FAIL                                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                            # Where to send mail	

# Load necessary modules (if jq is needed)
module load jq

# Set your input directory containing JSON files
input_dir="/scratch/ma95362/PRJNA823537_ET125/Tbprofiler_module/results"
output_dir="/scratch/ma95362/PRJNA823537_ET125/Tbprofiler_module/Parsing"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop over each JSON file in the directory
for json_file in "$input_dir"/*.json; do
    # Display the name of the current file being processed
    echo "Processing $json_file..."

    # Extract lineage information (e.g., major lineage and sublineage)
    lineage=$(jq -r '.lineage' "$json_file")
    sublineage=$(jq -r '.sublineage // "Not available"' "$json_file")

    # Extract spoligotype information
    spoligotype=$(jq -r '.spoligotype // "Not available"' "$json_file")

    # Extract drug resistance information for key drugs
    drug_resistance_isoniazid=$(jq -r '.drug_resistance.inh // "Not available"' "$json_file")
    drug_resistance_rifampicin=$(jq -r '.drug_resistance.rif // "Not available"' "$json_file")
    drug_resistance_streptomycin=$(jq -r '.drug_resistance.str // "Not available"' "$json_file")
    drug_resistance_ethambutol=$(jq -r '.drug_resistance.emb // "Not available"' "$json_file")
    drug_resistance_pza=$(jq -r '.drug_resistance.pza // "Not available"' "$json_file")
    drug_resistance_ethionamide=$(jq -r '.drug_resistance.eth // "Not available"' "$json_file")
    drug_resistance_capreomycin=$(jq -r '.drug_resistance.cap // "Not available"' "$json_file")
    drug_resistance_ofloxacin=$(jq -r '.drug_resistance.ofx // "Not available"' "$json_file")
    drug_resistance_amikacin=$(jq -r '.drug_resistance.amk // "Not available"' "$json_file")
    drug_resistance_kanamycin=$(jq -r '.drug_resistance.kan // "Not available"' "$json_file")
    
    # Prepare the output file name
    output_file="$output_dir/$(basename "$json_file" .json)_parsed.txt"

    # Save the extracted information to a new text file
    echo "Lineage: $lineage" > "$output_file"
    echo "Sublineage: $sublineage" >> "$output_file"
    echo "Spoligotype: $spoligotype" >> "$output_file"
    echo "Isoniazid Resistance: $drug_resistance_isoniazid" >> "$output_file"
    echo "Rifampicin Resistance: $drug_resistance_rifampicin" >> "$output_file"
    echo "Streptomycin Resistance: $drug_resistance_streptomycin" >> "$output_file"
    echo "Ethambutol Resistance: $drug_resistance_ethambutol" >> "$output_file"
    echo "PZA Resistance: $drug_resistance_pza" >> "$output_file"
    echo "Ethionamide Resistance: $drug_resistance_ethionamide" >> "$output_file"
    echo "Capreomycin Resistance: $drug_resistance_capreomycin" >> "$output_file"
    echo "Ofloxacin Resistance: $drug_resistance_ofloxacin" >> "$output_file"
    echo "Amikacin Resistance: $drug_resistance_amikacin" >> "$output_file"
    echo "Kanamycin Resistance: $drug_resistance_kanamycin" >> "$output_file"

    # Optional: Display the contents of the output file (for debugging purposes)
    cat "$output_file"
done

echo "Parsing completed!"
