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
PANAROO_DIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_local_all/panaroo"
TREEFILE=/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_local_all/iqtree/core-genome.contree          
PYSEER_OUT=/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_local_all/pyseer_output
METADATA="/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_local_all/metadata_local.tab"
PRES="/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_local_all/panaroo/filtered_output/gene_presence_absence_filt_pseudo_length.Rtab"

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
echo "Running Pyseer GWAS for EPTB vs Pulmonary"

    pyseer \
        --lmm \
        --phenotypes $METADATA \
        --phenotype-column $Phenotype_numeric \
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