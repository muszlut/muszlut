#!/bin/bash
#SBATCH --job-name=TBprofiler_Mamba         
#SBATCH --partition=batch                    
#SBATCH --ntasks=1                            
#SBATCH --cpus-per-task=32                     
#SBATCH --mem=120gb                            
#SBATCH --time=05-00:00:00                    
#SBATCH --output=/scratch/ma95362/all_in_all_reads/bactopia_prepare/my_tbprofiler_results/logs/log.%j.out   
#SBATCH --error=/scratch/ma95362/all_in_all_reads/bactopia_prepare/my_tbprofiler_results/logs/log.%j.err    
#SBATCH --mail-type=END,FAIL                  
#SBATCH --mail-user=ma95362@uga.edu          

set -euo pipefail  # safer bash settings

# Load Micromamba module
module load Micromamba/2.3.0

# Initialize Micromamba for bash
eval "$(micromamba shell hook --shell bash)"

# Activate TB-Profiler environment
micromamba activate tbprofiler

# Set working directories
FASTQ_DIR="/scratch/ma95362/all_in_all_reads"
OUTDIR="/scratch/ma95362/all_in_all_reads/bactopia_prepare/my_tbprofiler_results"
FOFN="/scratch/ma95362/all_in_all_reads/bactopia_prepare/samples.fofn"

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"
cd "$OUTDIR"
# Loop through each line in the FOFN (sample, read1, read2)
while read -r sample read1 read2; do
    echo "Processing $sample..."
    SAMPLE_OUT="$OUTDIR/$sample"
    mkdir -p "$SAMPLE_OUT"

    tb-profiler profile \
        -1 "$read1" \
        -2 "$read2" \
        -p "$sample" \
        --db /home/ma95362/.conda/envs/mtbvartools/share/tbprofiler/tbdb \
        --spoligotype \
        --dir "$SAMPLE_OUT" \
        --threads 32 \
        --txt
done < "$FOFN"

# Deactivate the Micromamba environment
micromamba deactivate