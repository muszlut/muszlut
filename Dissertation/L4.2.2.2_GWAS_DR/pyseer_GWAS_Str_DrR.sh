#!/bin/bash
#SBATCH --job-name=pyseer_gwas
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# ------------------------------
# 1. Activate Conda environment
# ------------------------------
module load Miniforge3
source activate pyseer-env

# ------------------------------
# 2. Define directories and files
# ------------------------------
PANAROO_DIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/pangenome_tools_results/bactopia/bactopia-runs/pangenome_of_L4.2.2.2.2/panaroo/filtered_output"
TREEFILE="/scratch/ma95362/eth_national_analysis/all_fastq_reads/pangenome_tools_results/bactopia/bactopia-runs/pangenome_of_L4.2.2.2.2/iqtree/core-genome.treefile"
PYSEER_OUT="${PANAROO_DIR}/pyseer_out"
METADATA="${PANAROO_DIR}/metadata_all_numeric.tab"
PRES="${PANAROO_DIR}/struct_presence_absence.Rtab"

# Ensure output directory exists
mkdir -p "${PYSEER_OUT}"
cd "${PANAROO_DIR}" || { echo "❌ ERROR: Panaroo directory not found."; exit 1; }

# ------------------------------
# 3. (Skipped) Phylogenetic distance generation
# ------------------------------
# You already generated phylogeny_K.tsv manually

# ------------------------------
# 4. Run GWAS analysis
# ------------------------------
echo "Running Pyseer GWAS for Structural variant analysis for DR..."

# Define the phenotype column (must match column name in metadata_local.tab)
PHENOCOL="Dr_type_binary"

# Run Pyseer LMM
pyseer \
    --lmm \
    --phenotypes "${METADATA}" \
    --phenotype-column "${PHENOCOL}" \
    --pres "${PRES}" \
    --similarity "${PYSEER_OUT}/phylogeny_K.tsv" \
    --cpu 16 \
    --output-patterns "${PYSEER_OUT}/STR_DR_patterns_DR.txt" \
    > "${PYSEER_OUT}/Str_DR_gwas.txt"

echo "✅ Pyseer GWAS completed for DR."

# ------------------------------
# 5. Completion message
# ------------------------------
echo "✅ All Pyseer analyses completed successfully on $(date)."

# ------------------------------
# 6. Deactivate environment
# ------------------------------
conda deactivate
