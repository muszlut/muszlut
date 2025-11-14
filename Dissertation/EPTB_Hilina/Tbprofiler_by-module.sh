#!/bin/bash
#SBATCH --job-name=TBprofiler_module
#SBATCH --partition=batch           
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32             
#SBATCH --mem=120G                    
#SBATCH --time=03-00:00:00
#SBATCH --output=/scratch/ma95362/EPTB_Hilina/Bactopia_Run/split_reads/TBprofiler_by_module/log.%j.out   
#SBATCH --error=/scratch/ma95362/EPTB_Hilina/Bactopia_Run/split_reads/TBprofiler_by_module/log.%j.err    
#SBATCH --mail-type=END,FAIL                 
#SBATCH --mail-user=ma95362@uga.edu         

module load TB-Profiler/6.6.5

# Set working directories
FASTQ_DIR="/scratch/ma95362/EPTB_Hilina/Bactopia_Run/split_reads"
OUTDIR="/scratch/ma95362/EPTB_Hilina/Bactopia_Run/split_reads/TBprofiler_by_module"
FOFN="/scratch/ma95362/EPTB_Hilina/Bactopia_Run/TBprofiler_reads.fofn"
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
        --threads 8 \
        --txt
done < "$FOFN"
