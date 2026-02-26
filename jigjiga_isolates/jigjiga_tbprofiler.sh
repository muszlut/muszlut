#!/bin/bash
#SBATCH --job-name=TBprofiler_Mamba
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --time=07-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

set -euo pipefail

module load Micromamba/2.3.0
eval "$(micromamba shell hook --shell bash)"
micromamba activate tbprofiler_env

READS_DIR="/scratch/ma95362/clean_sequences_reads"
OUTDIR="/scratch/ma95362/publication/tbprofiler_results"

mkdir -p "$OUTDIR"
cd "$READS_DIR"

for R1 in *_R1.fastq.gz
do
    SAMPLE=$(basename "$R1" _R1.fastq.gz)
    R2="${SAMPLE}_R2.fastq.gz"

    echo "Processing $SAMPLE"

    tb-profiler profile \
        --read1 "$R1" \
        --read2 "$R2" \
        --db /home/ma95362/.conda/envs/tbprofiler_env/share/tbprofiler/tbdb\
        --threads 16 \
        --spoligotype \
        --prefix "$SAMPLE" \
        --txt \
        --dir "$OUTDIR"
done