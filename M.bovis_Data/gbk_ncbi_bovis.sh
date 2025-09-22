#!/bin/bash
#SBATCH --job-name=M.bovis_ref_fasta_ncbi       # Job name
#SBATCH --partition=batch                        # Partition (queue) name
#SBATCH --ntasks=1                               # Run on a single CPU
#SBATCH --cpus-per-task=8                        # Number of cores per task
#SBATCH --mem=40gb                               # Job memory request
#SBATCH --time=03:00:00                          # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out   # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err    # Standard error log

#SBATCH --mail-type=END,FAIL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu              # Where to send mail

# Set output directory variable
OUTDIR="/scratch/ma95362/ETH_bovis_Sequence/bovis_REF/Fasta"

# Set URL for AF2122/97 reference genome (GCF_000195835.3)
URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/835/GCF_000195835.3_ASM19583v2/GCF_000195835.3_ASM19583v2_genomic.fna.gz"

# Create output directory if it doesn't exist
mkdir -p $OUTDIR
cd $OUTDIR

# Log start time
echo "Download started at $(date)"

# Download and decompress the reference genome
curl -s --fail $URL | gunzip -c > $OUTDIR/AF2122_97.fasta

# Check if download was successful
if [ $? -eq 0 ]; then
    echo "Download completed successfully at $(date)"
else
    echo "Download failed at $(date)" >&2
    exit 1
fi

# Optional: Preview the header
echo "Reference header:"
head -n 1 $OUTDIR/AF2122_97.fasta