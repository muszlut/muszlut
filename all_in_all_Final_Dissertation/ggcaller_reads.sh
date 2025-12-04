#!/bin/bash
#SBATCH --job-name=ggcaller_all_samples
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/ma95362/all_in_all_reads/logs/ggcaller_all.out
#SBATCH --error=/scratch/ma95362/all_in_all_reads/logs/ggcaller_all.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu
# Load environment
module purge
module load ggCaller/1.4.1
# Define paths
READS_LIST=/scratch/ma95362/all_in_all_reads/input_reads.txt   # list of all paired-end samples
REFS_LIST=/scratch/ma95362/all_in_all_reads/refs_list.txt    # curated reference genomes
BALROG_DB=/scratch/ma95362/ggcaller_db
OUTDIR=/scratch/ma95362/all_in_all_reads/ggcaller_results
# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"
# Move to working directory
cd $OUTDIR
# Run ggCaller on all 1500 samples together
ggcaller \
    --reads ${READS_LIST} \                  # list file with all paired-end reads
    --refs ${REFS_LIST} \                    # reference list for annotation consistency
    --annotation ultrasensitive \            # sensitive annotation mode
    --alignment pan \                        # pangenome alignment
    --aligner def \                          # default aligner
    --balrog-db ${BALROG_DB} \               # Balrog DB for clonal organisms
    --threads 24 \                           # use all allocated CPUs
    --save \                                 # save intermediate outputs
    --out ${OUTDIR}                          # unified output directory