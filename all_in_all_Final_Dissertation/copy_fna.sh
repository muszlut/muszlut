#!/bin/bash
#SBATCH --job-name=copy_fna
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/ma95362/ggcaller_module/logs/log.%j.out
#SBATCH --error=/scratch/ma95362/ggcaller_module/logs/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Paths
src="/scratch/ma95362/all_in_all_reads/bactopia_prepare"
dst="/scratch/ma95362/all_in_all_reads/ggcaller_pangenome_fna"
exclude="${src}/bactopia-exclude.tsv"

mkdir -p "$dst"

echo "Reading exclude list..."
declare -A exclude_map

# Load exclude list into a hash map for fast lookup
while read -r sample; do
    exclude_map["$sample"]=1
done < "$exclude"

echo "Starting copy process..."

# Loop through all sample directories
for sample_dir in "$src"/*/; do
    sample=$(basename "$sample_dir")

    # Skip if sample is in exclude list
    if [[ ${exclude_map[$sample]} ]]; then
        echo "Skipping excluded sample: $sample"
        continue
    fi

    # Path to fna.gz
    fna="${sample_dir}/main/assembler/${sample}.fna.gz"

    # Copy only if the file exists
    if [[ -f "$fna" ]]; then
        echo "Copying $sample..."
        cp "$fna" "$dst/"
    else
        echo "Missing fna.gz for $sample, skipping."
    fi
done

echo "Copy completed."
