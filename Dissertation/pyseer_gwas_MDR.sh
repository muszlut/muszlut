#!/bin/bash
#SBATCH --job-name=pyseer_gwas_MDR
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
PYSEER_OUT="${PANAROO_DIR}/pyseer_drugtype_out"
METADATA="${PANAROO_DIR}/metadata_resistance_level.tab"
PRES="${PANAROO_DIR}/gene_presence_absence_filt_pseudo_length.Rtab"
K_MATRIX="${PANAROO_DIR}/pyseer_out/phylogeny_K.tsv"

# Create output directory
mkdir -p $PYSEER_OUT
cd $PANAROO_DIR || exit 1

# ------------------------------
# 3. Run GWAS
# ------------------------------
echo "Running Pyseer GWAS for drugtype"

pyseer \
    --lmm \
    --phenotypes $METADATA \
    --phenotype-column dr_type \
    --pres $PRES \
    --similarity $K_MATRIX \
    --cpu 16 \
    --output-patterns ${PYSEER_OUT}/resistance_level_gwas_patterns.txt \
    > ${PYSEER_OUT}/resistance_level_gwas.txt

echo "✅ Pyseer GWAS for drugtype completed successfully on $(date)"
conda deactivate

#To get the threshold and patterns, run the following command in the directory of pyseer output after the above script is done
#python ~/pyseer/scripts/count_patterns.py ./gene_patterns_drugtype.txt (match this the output pattern file name)
#Patterns:       412
#Threshold:      1.21E-04
