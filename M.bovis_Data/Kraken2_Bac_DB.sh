#!/bin/bash
#SBATCH --job-name=kraken2_db_build
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120GB
#SBATCH --time=07-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Load Kraken2 module
module load Kraken2/2.1.3-gompi-2023a

# Define database location
DB_DIR=/scratch/ma95362/kraken2_db

# Create DB directory if it does not exist
mkdir -p $DB_DIR
cd $DB_DIR || { echo "Cannot cd to $DB_DIR"; exit 1; }

# Step 1: Download taxonomy
kraken2-build --download-taxonomy --db $DB_DIR

# Step 2: Download bacterial genomes
kraken2-build --download-library bacteria --db $DB_DIR

# Step 3: Build the database
kraken2-build --build --db $DB_DIR --threads 16
