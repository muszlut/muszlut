#!/bin/bash
#SBATCH --job-name=Tbprofiler_FOFN
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120gb
#SBATCH --time=03-00:00:00
#SBATCH --output=/scratch/ma95362/EPTB_Hilina/Newe/logs/log.%j.out
#SBATCH --error=/scratch/ma95362/EPTB_Hilina/Newe/logs/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Run as an executable file
chmod +x generate_fofn.sh

# Set the directory with your reads
READ_DIR="/scratch/ma95362/EPTB_Hilina/Bactopia_Run/split_reads"

# Output file
FOFN="/scratch/ma95362/EPTB_Hilina/Bactopia_Run/TBprofiler_reads.fofn"

# Generate FOFN
for r1 in "$READ_DIR"/*_R1.fastq.gz; do
  base=$(basename "$r1" | sed 's/_R1.fastq.gz//' | sed 's/\.$//')
  r2="${r1/_R1.fastq.gz/_R2.fastq.gz}"
  echo -e "${base}\t${r1}\t${r2}"
done > "$FOFN"

echo "FOFN written to $FOFN"