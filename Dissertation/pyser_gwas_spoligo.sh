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
PYSEER_OUT="${PANAROO_DIR}/pyseer_T3ETH_out"
METADATA="${PANAROO_DIR}/meta_numeric_T3ETH.tab"
PRES="${PANAROO_DIR}/gene_presence_absence_filt_pseudo_length.Rtab"
K_MATRIX="${PANAROO_DIR}/pyseer_out/phylogeny_K.tsv"

# Create output directory
mkdir -p $PYSEER_OUT
cd $PANAROO_DIR || exit 1

# ------------------------------
# 3. Run GWAS
# ------------------------------
echo "Running Pyseer GWAS for T3-ETH family vs others..."

pyseer \
    --lmm \
    --phenotypes $METADATA \
    --phenotype-column spo_family \
    --pres $PRES \
    --similarity $K_MATRIX \
    --cpu 16 \
    --output-patterns ${PYSEER_OUT}/gene_patterns_T3ETH.txt \
    > ${PYSEER_OUT}/T3ETH_gwas.txt

echo "✅ Pyseer GWAS for T3-ETH completed successfully on $(date)"
conda deactivate