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
METADATA="${PANAROO_DIR}/metadata.tab"
PRES="${PANAROO_DIR}/gene_presence_absence_filt_pseudo_length.Rtab"

# Create output directory
mkdir -p $PYSEER_OUT
cd $PANAROO_DIR || exit 1

# ------------------------------
# 3. Generate phylogenetic distance matrix
# ------------------------------
echo "Running phylogenetic distance matrix generation..."
phylogeny_distance.py --lmm $TREEFILE > ${PYSEER_OUT}/phylogeny_K.tsv

# ------------------------------
# 4. Run GWAS for each antibiotic
# ------------------------------
echo "Starting GWAS for all antibiotics..."
for anti in rifampicin isoniazid ethambutol pyrazinamide moxifloxacin levofloxacin bedaquiline delamanid pretomanid linezolid streptomycin amikacin kanamycin capreomycin clofazimine ethionamide para-aminosalicylic_acid cycloserine
do
    echo "Running Pyseer for ${anti}..."
    pyseer \
        --lmm \
        --phenotypes $METADATA \
        --phenotype-column $anti \
        --pres $PRES \
        --similarity ${PYSEER_OUT}/phylogeny_K.tsv \
        --cpu 16 \
        --output-patterns ${PYSEER_OUT}/gene_patterns_${anti}.txt \
        > ${PYSEER_OUT}/${anti}_gwas.txt
done

# ------------------------------
# 5. Completion message
# ------------------------------
echo "✅ Pyseer analysis completed successfully on $(date)"
# Deactivate Conda environment
conda deactivate