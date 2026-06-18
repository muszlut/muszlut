#!/bin/bash
#SBATCH --job-name=bakta_db
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=80G
#SBATCH --time=07-00:00:00
#SBATCH --output=/scratch/ma95362/bakta_logs/bakta_db_%j.out
#SBATCH --error=/scratch/ma95362/bakta_logs/bakta_db_%j.err

mkdir -p /scratch/ma95362/bakta_logs
mkdir -p /scratch/ma95362/full_bakta_db
cd /scratch/ma95362/full_bakta_db

wget -c https://zenodo.org/record/14916843/files/db.tar.xz

xz -t db.tar.xz

tar -xJf db.tar.xz