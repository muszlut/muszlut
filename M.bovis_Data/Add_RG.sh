#!/bin/bash
#SBATCH --job-name=add_readgroups
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/ma95362/logs/add_readgroups_%j.out
#SBATCH --error=/scratch/ma95362/logs/add_readgroups_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Load conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bovisanalyzer

# Directory containing your work subfolders
WORKDIR="/scratch/ma95362/ETH_bovis_Sequence/bovisanalyzer_output/2nd_run/work"

# Loop through all sorted BAMs
find $WORKDIR -name "*.sorted.bam" | while read BAM; do
    FILENAME=$(basename "$BAM" .sorted.bam)
    DIRNAME=$(dirname "$BAM")
    OUTBAM="${DIRNAME}/${FILENAME}.rg.bam"

    echo "Processing $BAM ..."

    picard AddOrReplaceReadGroups \
        I="$BAM" \
        O="$OUTBAM" \
        RGID="$FILENAME" \
        RGLB=lib1 \
        RGPL=ILLUMINA \
        RGPU=unit1 \
        RGSM="$FILENAME" \
        CREATE_INDEX=true

    echo "Done: $OUTBAM"
done

# Deactivate conda
conda deactivate
