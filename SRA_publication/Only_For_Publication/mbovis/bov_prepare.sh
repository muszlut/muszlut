#!/bin/bash
#SBATCH --job-name=Bactopia_bovis
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/ma95362/Bio_project_publication/Mbovis/logs/log.%j.out
#SBATCH --error=/scratch/ma95362/Bio_project_publication/Mbovis/logs/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

set -euo pipefail

# Load Mamba
module load Mamba/23.11.0-0

# Activate Bactopia environment
source activate bactopia4

# -------------------------------
# Variables
# -------------------------------
READS_DIR="/scratch/ma95362/Bio_project_publicatio/mbovis_seq"
OUTDIR="/scratch/ma95362/Bio_project_publication/Mbovis"
#REF="/scratch/ma95362/gbk/ncbi_dataset/data/GCF_000195955.2/genomic.gbk"

# -------------------------------
# Create output directory
# -------------------------------
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# -------------------------------
# Prepare samples
# -------------------------------
#bactopia prepare \
#    --path "$READS_DIR" \
#    --species "Mycobacterium bovis" \
#    --genome-size 4345492 \
#    > "$OUTDIR/samples.txt"
## -------------------------------
# Run Bactopia
# -------------------------------
export NXF_OPTS="-Xms1g -Xmx6g"
bactopia \
    --samples "$OUTDIR/samples.txt" \
    --coverage 100 \
    --outdir "$OUTDIR" \
    -c sapelo2.config \
    -work-dir /scratch/ma95362/Bio_project_publication/Mbovis/work \
    -resume
bactopia summary --bactopia-path "$OUTDIR"
# Run Snippy
# Nextflow JVM memory
#export NXF_OPTS="-Xms2g -Xmx8g"
#Then, now run SNippy workflow with the reference genome and the Bactopia output directory
#bactopia \
#    --wf snippy \
#    --reference "$REF" \
#    --bactopia "$OUTDIR" \
#    --exclude snippy_exclude.txt \
#    --outdir /scratch/ma95362/Bio_project_publication/snippy_results \
#    --max_cpus 24 \
#    -process.maxForks 8 \
#    -resume
#``
#bactopia \
#    --wf pangenome \
#    --bactopia $OUTDIR