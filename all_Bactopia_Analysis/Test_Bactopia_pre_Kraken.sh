#!/bin/bash
#SBATCH --job-name=M.bovis_Bactopia                        # Job name
#SBATCH --partition=batch                                   # Partition
#SBATCH --ntasks=1                                          # Single task
#SBATCH --cpus-per-task=8                                   # CPUs per task
#SBATCH --mem=40gb                                         # Memory
#SBATCH --time=05-00:00:00                                  # Walltime
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out       # STDOUT
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err        # STDERR
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

#----------------------------
# Directories
#----------------------------
OUTDIR="/scratch/ma95362/ETH_bovis_Sequence/Bactopia_Analysis"
mkdir -p $OUTDIR
cd $OUTDIR

#----------------------------
# Modules
#----------------------------
module load Bactopia/3.2.0

# Prevent Python from using user-local packages
export PYTHONNOUSERSITE=1

# Activate Bactopia conda environment and ensure compatible NumPy/SciPy
source /home/ma95362/.bactopia/conda/bioconda--bactopia-sketcher-1.0.1/bin/activate
conda install -y -c conda-forge numpy=1.26 scipy

#----------------------------
# Prepare samples
#----------------------------
bactopia prepare \
    --path /scratch/ma95362/ETH_bovis_Sequence/reads/64_SH_Sequence_data/raw/sub_raw \
    --species "Mycobacterium tuberculosis" \
    --genome-size 4410000 \
    > $OUTDIR/MGA_samples.txt

#----------------------------
# Run Bactopia
#----------------------------
bactopia \
    --samples $OUTDIR/MGA_samples.txt \
    --coverage 100 \
    --outdir $OUTDIR/MGA_paired_end_samples \
    --max_cpus 8 \
    --skip_qc false \
    --skip_kraken2 false

#----------------------------
# Generate summary
#----------------------------
bactopia summary \
    --bactopia-path $OUTDIR/MGA_paired_end_samples \
    --outdir $OUTDIR/MGA_summary

#----------------------------
# Extract top 5 Kraken2 hits and export to Excel
#----------------------------
python3 <<EOF
import pandas as pd
import glob
import os

kraken_files = glob.glob("$OUTDIR/MGA_paired_end_samples/*/kraken2/*_report.txt")
all_data = []

for f in kraken_files:
    sample_name = os.path.basename(os.path.dirname(os.path.dirname(f)))
    df = pd.read_csv(f, sep="\t", header=None, names=["perc", "reads", "reads_clade", "rank", "taxid", "name"])
    df = df.sort_values(by="perc", ascending=False).head(5)  # top 5
    df["sample"] = sample_name
    all_data.append(df)

if all_data:
    final_df = pd.concat(all_data)
    final_df = final_df[["sample", "perc", "reads", "reads_clade", "rank", "taxid", "name"]]
    final_df.to_excel(os.path.join("$OUTDIR", "Top5_Kraken2.xlsx"), index=False)
EOF

echo "Bactopia pipeline finished. Top5 Kraken2 results saved to Top5_Kraken2.xlsx"
