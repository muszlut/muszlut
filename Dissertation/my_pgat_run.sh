#!/bin/bash
#SBATCH --job-name=sharma_pgat
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

module load anaconda3
source activate sharma_pgat

# path to your Rtab file and output dir
RTAB="/scratch/ma95362/eth_national_analysis/all_fastq_reads/pangenome_tools_results/bactopia/bactopia-runs/pangenome-20251011-115307/panaroo/gene_presence_absence.Rtab"
OUTDIR="/scratch/ma95362/eth_national_analysis/all_fastq_reads/pangenome_tools_results/bactopia/bactopia-runs/pangenome-20251011-115307/panaroo/output_pgat"
# Create output folder if not exist
mkdir -p $OUTDIR
cd $OUTDIR
# run pangenome analysis tool
python /home/ma95362/PanGenomeAnalysisTool/pan_genome_analysis.py -f $RTAB -i 100 -o $OUTDIR
source deactivate