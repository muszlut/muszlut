#!/bin/bash
#SBATCH --job-name=gunzip_job          # Job name
#SBATCH --partition=batch
#SBATCH --ntasks=1                     # One task
#SBATCH --cpus-per-task=4              # Number of CPU cores
#SBATCH --mem=8G                       # Memory (adjust as needed)
#SBATCH --time=02:00:00                # Max runtime (hh:mm:ss)
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# -----------------------------
# Load any required modules
# -----------------------------
#module load gzip/1.12   # or module load gzip (if available)
cd /scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_of_1368/panaroo
# -----------------------------
# Define input and output files
# -----------------------------
# Example 1: decompress a single file
gunzip -c gene_data.csv.gz > gene_data.csv

# Example to group all the significant hits and put it into a new csv file to show you which isolates have them
#grep -E "group_2728|group_2726|group_2586|group_2104" gene_presence_absence.csv > selected_high_significant_genes.csv


