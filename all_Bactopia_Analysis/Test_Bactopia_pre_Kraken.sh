#!/bin/bash
#SBATCH --job-name=Kraken2_Test_BacPrep                          # Job name
#SBATCH --partition=batch                                        # Partition (queue) name
#SBATCH --ntasks=1                                               # Run on a single CPU
#SBATCH --cpus-per-task=8                                        # Number of cores per task
#SBATCH --mem=40gb                                               # Job memory request
#SBATCH --time=05-00:00:00                                       # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out             # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err              # Standard error log

#SBATCH --mail-type=END,FAIL                                     # Mail events
#SBATCH --mail-user=ma95362@uga.edu                              # Where to send mail	

#----------------------------
# Set output directory
#----------------------------
OUTDIR="/scratch/ma95362/musse_MGA/all_reads_Bactopia_Analysis/Test_Bactopia_pre_Kraken"

# Make sure output directory exists
mkdir -p $OUTDIR

#----------------------------
# Load module
#----------------------------
module load Bactopia/3.2.0
cd $OUTDIR

#----------------------------
# Prepare sample sheet
#----------------------------
bactopia prepare \
    --path /scratch/ma95362/ETH_bovis_Sequence/reads/64_SH_Sequence_data/raw/sub_raw \
    --species "Mycobacterium tuberculosis variety bovis" \
    --genome_size 4200000-4500000 \
    > $OUTDIR/MGA_samples.txt

#----------------------------
# Run Bactopia pipeline
#----------------------------
bactopia \
    --samples $OUTDIR/MGA_samples.txt \
    --coverage 100 \
    --outdir $OUTDIR/MGA_paired_end_samples \
    --max_cpus 8 \
    --skip_qc false \
    --skip_kraken2 false

#----------------------------
# Generate summary report
#----------------------------
bactopia summary \
    --bactopia-path $OUTDIR/MGA_paired_end_samples

#----------------------------
# Collect top 5 Kraken2 hits
#----------------------------
TOPKRAKEN="$OUTDIR/top5_kraken2_summary.tsv"
echo -e "Sample\tRank\tPercent\tReads\tTaxonomy_ID\tRank_Code\tScientific_Name" > $TOPKRAKEN

for KREPORT in $OUTDIR/MGA_paired_end_samples/*/kraken2/*.report.txt; do
    SAMPLE=$(basename $(dirname $(dirname $KREPORT)))
    # Grab the first 5 lines, add a numeric rank (1–5)
    awk -v SAMP="$SAMPLE" 'NR<=5 {print SAMP"\t"NR"\t"$0}' OFS="\t" $KREPORT >> $TOPKRAKEN
done

#----------------------------
# Convert TSV to Excel
#----------------------------
module load Python/3.9.6-GCCcore-11.2.0
python3 - <<'EOF'
import pandas as pd
import os
outdir = os.environ["OUTDIR"]
tsv_file = os.path.join(outdir, "top5_kraken2_summary.tsv")
excel_file = os.path.join(outdir, "top5_kraken2_summary.xlsx")
df = pd.read_csv(tsv_file, sep="\t")
df.to_excel(excel_file, index=False)
EOF
#----------------------------

