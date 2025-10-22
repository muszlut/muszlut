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

# Go to your Panaroo filtered directory
cd /scratch/ma95362/eth_national_analysis/all_fastq_reads/pangenome_tools_results/bactopia/bactopia-runs/pangenome_of_L4.2.2.2.2/panaroo/filtered_output

# Create output directory
mkdir -p pyseer_out

# Step 1: Calculate phylogenetic distances
phylogeny_distance.py --lmm /scratch/ma95362/eth_national_analysis/all_fastq_reads/pangenome_tools_results/bactopia/bactopia-runs/pangenome_of_L4.2.2.2.2/iqtree/core-genome.treefile > pyseer_out/phylogeny_K.tsv

# Step 2: Run GWAS for each antibiotic phenotype
for anti in rifampicin isoniazid ethambutol pyrazinamide moxifloxacin levofloxacin bedaquiline delamanid pretomanid linezolid streptomycin amikacin kanamycin capreomycin clofazimine ethionamide para-aminosalicylic_acid cycloserine
do
  pyseer-runner.py --lmm \
  --phenotypes ./metadata.tab \
  --pres ./gene_presence_absence_filt_pseudo_length.Rtab \
  --similarity ./pyseer_out/phylogeny_K.tsv \
  --phenotype-column $anti \
  --output-patterns ./pyseer_out/gene_patterns_${anti}.txt \
  > ./pyseer_out/${anti}_gwas.txt
done

echo "✅ Pyseer GWAS completed successfully."
conda deactivate
