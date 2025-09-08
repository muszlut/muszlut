#!/bin/bash
#SBATCH --job-name=kraken2_merge
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=40gb
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/ma95362/scratch/merge.%j.out
#SBATCH --error=/scratch/ma95362/scratch/merge.%j.err

module load Python/3.11.5
KRAKENTOOLS=/scratch/ma95362/KrakenTools

WORKDIR=/scratch/ma95362/ETH_bovis_Sequence/reads/64_SH_Sequence_data/raw/fixed_fastqs/fixed_fastqs
cd $WORKDIR

# Collect reports and sample names in sorted order
REPORTS=$(ls *_kraken2_report.txt | sort)
SAMPLES=$(ls *_R1.fastq.gz | sort | sed 's/_R1.fastq.gz//')

# Run KrakenTools combine script
python $KRAKENTOOLS/combine_kreports.py \
    -r $REPORTS \
    -o kraken2_combined.tsv \
    --sample-names $SAMPLES \
    --display-headers
