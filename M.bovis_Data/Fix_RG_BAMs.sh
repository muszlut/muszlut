#!/bin/bash
#SBATCH --job-name=Fix_RG_BAMs
#SBATCH --partition=batch                                
#SBATCH --ntasks=1                                       
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Load conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bovisanalyzer


# Load required modules
module load picard/2.27.5
module load samtools

# Directory with BAMs (edit this if needed)
BAM_DIR="/scratch/ma95362/ETH_bovis_Sequence/bovisanalyzer_output/2nd_run/work"
cd "$BAM_DIR" || exit 1

echo "=== Starting Read Group Fixing Job ==="
echo "Working directory: $PWD"
echo "Start time: $(date)"
echo "--------------------------------------"

# Loop through all BAM files
for bam in *.bam; do
    # Skip if file does not exist
    [ -e "$bam" ] || continue

    # Check if BAM already has read group (RG) info
    if ! samtools view -H "$bam" | grep -q "^@RG"; then
        sample=$(basename "$bam" .bam)
        fixed="${sample}_RGfixed.bam"

        echo "🧩 Fixing missing read group for: $sample"

        # Add a standard read group
        picard AddOrReplaceReadGroups \
            I="$bam" \
            O="$fixed" \
            RGID="$sample" \
            RGLB="lib1" \
            RGPL="ILLUMINA" \
            RGPU="unit1" \
            RGSM="$sample" \
            VALIDATION_STRINGENCY=SILENT

        # Replace old BAM if successful
        if [ -s "$fixed" ]; then
            mv "$fixed" "$bam"
            echo "✔ Patched: $bam"
        else
            echo "❌ Failed to fix: $bam"
        fi
    else
        echo "✅ Already has read group: $bam"
    fi
done

echo "--------------------------------------"
echo "All BAMs processed. End time: $(date)"
echo "======================================"

# Deactivate conda environment
conda deactivate
echo "Conda environment deactivated."