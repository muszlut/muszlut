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
src="/scratch/ma95362/all_in_all_reads/ggcaller_pangenome_fna"
dst="/scratch/ma95362/all_in_all_reads/ggcaller_pangenome_fna_somali_only"

mkdir -p "$dst"
cd "$src"
# Loop through all .fna.gz files in source
for fna in "$src"/*.fna.gz; do
    sample=$(basename "$fna" .fna.gz)

    # Skip files starting with ERR or SRR
    if [[ "$sample" == ERR* || "$sample" == SRR* ]]; then
        echo "Skipping $sample..."
        continue
    fi

    # Copy only if the file exists
    if [[ -f "$fna" ]]; then
        echo "Copying $sample..."
        cp "$fna" "$dst/"
    else
        echo "Missing fna.gz for $sample, skipping."
    fi
done