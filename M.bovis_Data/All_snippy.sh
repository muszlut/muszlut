#!/bin/bash
#SBATCH --job-name=BovAll_Snippy
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# -----------------------------
# Load environment modules
# -----------------------------
module load Bactopia/3.2.0-conda
module load Java/17.0.6


#Set output directory variable
OUTDIR="/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads"
REF="/scratch/ma95362/ETH_bovis_Sequence/bovis_REF/genomic.gbk"

# Create output folder if not exist
mkdir -p $OUTDIR
cd $OUTDIR
# Run snippy workflow using all samples
bactopia \
    --wf snippy \
    --reference "$REF" \
    --bactopia "$OUTDIR/M.bovis_paired_end_samples" \
    --exclude /scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples/all_fastqs/bactopia-exclude.tsv