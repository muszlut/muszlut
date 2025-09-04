#!/bin/bash
#SBATCH --job-name=updated_Bactopia_Prep                 # Job name
#SBATCH --partition=batch                                # Partition (queue)
#SBATCH --ntasks=1                                       # Single task
#SBATCH --cpus-per-task=8                                # CPUs per task
#SBATCH --mem=40gb                                       # Memory
#SBATCH --time=05-00:00:00                               # Time limit (HH:MM:SS)
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out    # STDOUT log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err     # STDERR log
#SBATCH --mail-type=END,FAIL                             # Mail events
#SBATCH --mail-user=ma95362@uga.edu                      # Your email

#----------------------------
# Set output directory
#----------------------------
OUTDIR="/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis"
mkdir -p $OUTDIR

#----------------------------
# Load Bactopia module
#----------------------------
module load Bactopia/3.2.0
cd $OUTDIR

#----------------------------
# Prepare sample sheet
#----------------------------
bactopia prepare \
    --path /scratch/ma95362/ETH_bovis_Sequence/reads/64_SH_Sequence_data/raw \
    --species "Mycobacterium bovis" \
    --genome-size 4410000 \
    > $OUTDIR/M.bovis_samples.txt

#----------------------------
# Run Bactopia pipeline
#----------------------------
bactopia \
    --samples $OUTDIR/M.bovis_samples.txt \
    --coverage 100 \
    --outdir $OUTDIR/M.bovis_paired_end_samples \
    --max_cpus 8 \
    --skip_qc false \       # Ensure QC (Kraken2, FastQC, etc.) runs
    --skip_kraken2 false    # Explicitly run Kraken2 for species ID

#----------------------------
# Generate summary report
#----------------------------
bactopia summary \
    --bactopia-path $OUTDIR/M.bovis_paired_end_samples \
    --outdir $OUTDIR/M.bovis_summary

#----------------------------
# Optional: Extract top Kraken2 species hits for all samples
#----------------------------
echo -e "Sample\tTop_Kraken2_Species" > $OUTDIR/M.bovis_species_summary.txt
for f in $OUTDIR/M.bovis_paired_end_samples/*/tools/kraken2/*.report; do
    sample=$(basename $(dirname $(dirname $f)))
    top_species=$(awk '{print $6; exit}' $f)
    echo -e "${sample}\t${top_species}" >> $OUTDIR/M.bovis_species_summary.txt
done
