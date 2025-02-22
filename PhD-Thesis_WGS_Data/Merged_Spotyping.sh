#!/bin/bash
#SBATCH --job-name=Musse_Merged_SpoTyping
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40gb
#SBATCH --time=05-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Set absolute paths (Slurm jobs often run in a temporary working directory)
OUTDIR="/scratch/ma95362/musse_MGA/merged/SpoTyping_results"
mkdir -p "${OUTDIR}"

# Load modules
module purge  # Start with a clean environment
module load Bactopia/3.1.0
module load Singularity/3.11.4-GCCcore-11.3.0
module load conda

# --- Step 1: Run Bactopia ---
# Use absolute paths for inputs/outputs
bactopia \
  --samples samplesheet.csv \
  --path /work/fdqlab/merged_reads/Ethiopia_wgs_mtb_2024/first_run \
  --outdir "${OUTDIR}/bactopia_results" \  # Changed to absolute path
  --wf myco \
  --genome_size 4.4M \
  --species "Mycobacterium tuberculosis" \
  -profile singularity

# Check if Bactopia succeeded
if [ ! -d "${OUTDIR}/bactopia_results/assemblies" ]; then
  echo "ERROR: Bactopia assemblies not found!"
  exit 1
fi

# --- Step 2: SpoTyping Setup ---
# Initialize conda and use absolute paths
source $(conda info --base)/etc/profile.d/conda.sh
conda create -n spotyping -y python=3.8 blast=2.13.0
conda activate spotyping

cd "${OUTDIR}"  # Work within the output directory
git clone https://github.com/xiaeryu/SpoTyping.git

# --- Step 3: Run SpoTyping ---
mkdir -p "${OUTDIR}/spotyping_results"

for ASSEMBLY in "${OUTDIR}/bactopia_results/assemblies/"*.fna.gz; do
  # Handle spaces in filenames (good practice)
  BASENAME=$(basename "${ASSEMBLY}" .fna.gz)
  
  # Decompress to a temporary directory
  mkdir -p "${OUTDIR}/temp"
  gunzip -c "${ASSEMBLY}" > "${OUTDIR}/temp/${BASENAME}.fna"
  
  # Run SpoTyping with full path
  python "${OUTDIR}/SpoTyping/SpoTyping.py" \
    -s "${OUTDIR}/temp/${BASENAME}.fna" \
    -o "${OUTDIR}/spotyping_results/${BASENAME}"
  
  # Cleanup
  rm -rf "${OUTDIR}/temp"
done

# --- Step 4: Aggregate Results ---
if [ -n "$(ls "${OUTDIR}/spotyping_results/"*"/spotyping.txt" 2>/dev/null)" ]; then
  echo "Sample,Spoligotype,SIT" > "${OUTDIR}/spotyping_results/spotyping_summary.csv"
  for RESULT in "${OUTDIR}/spotyping_results/"*"/spotyping.txt"; do
    SAMPLE=$(basename "$(dirname "${RESULT}")")
    SPOLIGOTYPE=$(awk '/Binary/ {print $3}' "${RESULT}")
    SIT=$(awk '/SIT/ {print $3}' "${RESULT}")
    echo "${SAMPLE},${SPOLIGOTYPE},${SIT}" >> "${OUTDIR}/spotyping_results/spotyping_summary.csv"
  done
else
  echo "ERROR: No spoligotype results found!"
  exit 1
fi