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

# Activate conda environment
source ~/.bashrc
conda activate pyseer-env

# Define paths
PANAROO_DIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/pangenome_tools_results/bactopia/bactopia-runs/pangenome_of_L4.2.2.2.2/panaroo/filtered_output"
TREEFILE="/scratch/ma95362/eth_national_analysis/all_fastq_reads/pangenome_tools_results/bactopia/bactopia-runs/pangenome_of_L4.2.2.2.2/iqtree/core-genome.treefile"
PYSEER_OUT="${PANAROO_DIR}/pyseer_out"
METADATA="${PANAROO_DIR}/metadata.tab"
PRES="${PANAROO_DIR}/gene_presence_absence_filt_pseudo_length.Rtab"

# Create output directory
mkdir -p $PYSEER_OUT

echo "Running phylogenetic distance matrix generation..."
python ~/pyseer/scripts/phylogeny_distance.py --lmm $TREEFILE > ${PYSEER_OUT}/phylogeny_K.tsv

# List of antibiotic phenotype columns (adjust based on your metadata)
for anti in rifampicin isoniazid ethambutol pyrazinamide moxifloxacin levofloxacin bedaquiline delamanid pretomanid linezolid streptomycin amikacin kanamycin capreomycin clofazimine ethionamide para-aminosalicylic_acid cycloserine
do
    echo "Running Pyseer for ${anti}..."
    python ~/pyseer/pyseer-runner.py \
        --lmm \
        --phenotypes $METADATA \
        --pres $PRES \
        --similarity ${PYSEER_OUT}/phylogeny_K.tsv \
        --phenotype-column $anti \
        --output-patterns ${PYSEER_OUT}/gene_patterns_${anti}.txt \
        > ${PYSEER_OUT}/${anti}_gwas.txt
done

echo "Pyseer analysis completed successfully!"
