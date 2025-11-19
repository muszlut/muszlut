#!/bin/bash
#SBATCH --job-name=pyseer_gwas
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/ma95362/eth_national_analysis/all_fastq_reads/pangenome_tools_results/bactopia/bactopia-runs/pangenome-20251117-164856/panaroo/panaroo_scratch_log/log.%j.out
#SBATCH --error=/scratch/ma95362/eth_national_analysis/all_fastq_reads/pangenome_tools_results/bactopia/bactopia-runs/pangenome-20251117-164856/panaroo/panaroo_scratch_log/log.%j.err
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
PANAROO_DIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/pangenome_tools_results/bactopia/bactopia-runs/pangenome-20251117-164856/panaroo/filtered_output"
TREEFILE="/scratch/ma95362/eth_national_analysis/all_fastq_reads/pangenome_tools_results/bactopia/bactopia-runs/pangenome-20251117-164856/iqtree/core-genome.treefile"
PYSEER_OUT="${PANAROO_DIR}/pyseer_out"
METADATA="${PANAROO_DIR}/New_L4_metadata.txt"
PRES="${PANAROO_DIR}/gene_presence_absence_filt_pseudo_length.Rtab"

# Create output directory
mkdir -p $PYSEER_OUT
cd $PANAROO_DIR || exit 1

# ------------------------------
# 3. Generate phylogenetic distance matrix
# ------------------------------
#echo "Running phylogenetic distance matrix generation..."
#phylogeny_distance.py --lmm $TREEFILE > ${PYSEER_OUT}/phylogeny_K.tsv
#I just did the above commands and generate phylogeny_K.tsv file manually then continue below
# ------------------------------
# 4. Run GWAS for each antibiotic
# ------------------------------
echo "Running Pyseer GWAS for DR..."

# Define the phenotype column (must match column name in metadata_local.tab)
PHENOCOL="Spoligo_T3_ETHfamily"

# Run Pyseer LMM
pyseer \
    --lmm \
    --phenotypes "${METADATA}" \
    --phenotype-column "${PHENOCOL}" \
    --pres "${PRES}" \
    --similarity "${PYSEER_OUT}/phylogeny_K.tsv" \
    --cpu 16 \
    --output-patterns "${PYSEER_OUT}/gene_patterns_New_L4_T3_ETHfamily.txt" \
    > "${PYSEER_OUT}/New_L4_T3_ETHfamily_gwas.txt"

echo "✅ Pyseer GWAS completed for New_L4_T3_ETHfamily."

# ------------------------------
# 5. Completion message
# ------------------------------
echo "✅ All Pyseer analyses completed successfully on $(date)."

# ------------------------------
# 6. Deactivate environment
# ------------------------------
conda deactivate