#!/bin/bash
#SBATCH --job-name=bovis_analyzer
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=05-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu


# Load conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bovisanalyzer

# Define inputs
SAMPLESHEET=/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples/samplesheet.csv
REFERENCE=/scratch/ma95362/ETH_bovis_Sequence/bovis_REF/Fasta/AF2122_97.fasta
KRAKEN2DB=/scratch/ma95362/kraken2_db/mini_db
OUTDIR=/scratch/ma95362/ETH_bovis_Sequence/bovisanalyzer_output

# Create output and logs directory
mkdir -p $OUTDIR
cd $OUTDIR
# add_read_groups.sh
# Add missing read groups to BAM files using Picard


for BAM in $(find /scratch/ma95362/ETH_bovis_Sequence/bovisanalyzer_output/work -name "*.sorted.bam"); do
    SAMPLE=$(basename $BAM .sorted.bam)
    OUTDIR=$(dirname $BAM)
    OUT="${OUTDIR}/${SAMPLE}.rg.bam"

    echo "Checking read groups for: $SAMPLE"

    if ! samtools view -H "$BAM" | grep -q "^@RG"; then
        echo "→ Adding read groups..."
        picard AddOrReplaceReadGroups \
            I="$BAM" \
            O="$OUT" \
            RGID="$SAMPLE" \
            RGLB="lib1" \
            RGPL="ILLUMINA" \
            RGPU="unit1" \
            RGSM="$SAMPLE" \
            CREATE_INDEX=true \
            VALIDATION_STRINGENCY=LENIENT
        mv "$OUT" "$BAM"
        echo "✔ Done for $SAMPLE"
    else
        echo "✔ Read groups already exist for $SAMPLE"
    fi
done
