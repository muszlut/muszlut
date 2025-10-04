#!/bin/bash
#SBATCH --job-name=bovis_analyzer
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=05-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

source $(conda info --base)/etc/profile.d/conda.sh
conda activate bovisanalyzer

SAMPLESHEET=/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples/samplesheet.csv
REFERENCE=/scratch/ma95362/ETH_bovis_Sequence/bovis_REF/Fasta/AF2122_97.fasta
KRAKEN2DB=/scratch/ma95362/mini_db        # <-- mini DB
OUTDIR=/scratch/ma95362/ETH_bovis_Sequence/bovisanalyzer_output

mkdir -p $OUTDIR
cd $OUTDIR

nextflow run avantonder/bovisanalyzer \
    -profile conda \
    --input $SAMPLESHEET \
    --reference $REFERENCE \
    --kraken2db $KRAKEN2DB \
    --outdir $OUTDIR \
    -resume \
    -c /scratch/ma95362/bovisanalyzer_custom.config


