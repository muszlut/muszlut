#!/bin/bash
#SBATCH --job-name=kraken2_db_build
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=250gb
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err            # Standard error log

#SBATCH --mail-type=END,FAIL                                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                            # Where to send mail	

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# -------------------------------
# Load module
# -------------------------------
module load Kraken2/2.1.5-gompi-2023a

# -------------------------------
# Set database directory
# -------------------------------
DBDIR="/scratch/ma95362/kraken2_db_clean"

mkdir -p "$DBDIR"
cd "$DBDIR"

# -------------------------------
# OPTIONAL: Prevent restart from scratch if interrupted
# -------------------------------
export KRAKEN2_DB="$DBDIR"

# -------------------------------
# Build STANDARD Kraken2 database
# -------------------------------
kraken2-build \
    --standard \
    --db "$DBDIR" \
    --threads 32 \
    --use-ftp