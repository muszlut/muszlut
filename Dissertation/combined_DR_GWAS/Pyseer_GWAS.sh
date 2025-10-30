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
PANAROO_DIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_of_1368/panaroo"
TREEFILE="/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_of_1368/core_tree.treefile"
PYSEER_OUT="/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_of_1368/panaroo/pyseer_DR_output"
METADATA="/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_of_1368/Full_metadata.tab"
PRES="/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_of_1368/panaroo/filtered_output/gene_presence_absence_filt_pseudo_length.Rtab"

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
echo "Running Pyseer GWAS for DR..."

# Define the phenotype column (must match column name in metadata_local.tab)
PHENOCOL="dr_classN"

# Run Pyseer LMM
pyseer \
    --lmm \
    --phenotypes "${METADATA}" \
    --phenotype-column "${PHENOCOL}" \
    --pres "${PRES}" \
    --similarity "${PYSEER_OUT}/phylogeny_K.tsv" \
    --cpu 16 \
    --output-patterns "${PYSEER_OUT}/gene_patterns_DR.txt" \
    > "${PYSEER_OUT}/DR_gwas.txt"

echo "✅ Pyseer GWAS completed for DR."

# ------------------------------
# 5. Completion message
# ------------------------------
echo "✅ All Pyseer analyses completed successfully on $(date)."

# ------------------------------
# 6. Deactivate environment
# ------------------------------
conda deactivate
