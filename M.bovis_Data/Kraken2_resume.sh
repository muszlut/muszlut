#!/bin/bash
#SBATCH --job-name=kraken2_resume
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=02-00:00:00
#SBATCH --output=/scratch/ma95362/kraken2_resume.%j.out
#SBATCH --error=/scratch/ma95362/kraken2_resume.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Load Conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bovisanalyzer

# Path
KRAKEN2DB=/scratch/ma95362/kraken2_db

# Continue/rebuild database.kraken
kraken2-build --build --db $KRAKEN2DB --threads 16 --fast-build

echo "database.kraken finished!"
```
