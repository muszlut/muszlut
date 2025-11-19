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
#cd /scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_of_1368/panaroo
cd /scratch/ma95362/eth_national_analysis/all_fastq_reads/pangenome_tools_results/bactopia/bactopia-runs/pangenome-20251117-164856/
# -----------------------------
# Define input and output files
# -----------------------------
# Example 1: decompress a single file
#gunzip -c gene_data.csv.gz > gene_data.csv 
#gunzip -c final_graph.gml.gz > final_graph.gml
#gunzip -c pan_genome_reference.fa.gz > pan_genome_reference.fa
#gunzip -c combined_DNA_CDS.fasta.gz > combined_DNA_CDS.fasta
#gunzip -c combined_protein_CDS.fasta.gz > combined_protein_CDS.fasta
gunzip -c core_gene_alignment.aln.gz > core_gene_alignment.aln
gunzip -c core-genome.aln.gz > core-genome.aln
gunzip -c core_gene_alignment_filtered.aln.gz > core_gene_alignment_filtered.aln