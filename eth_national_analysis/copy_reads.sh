#!/bin/bash
#SBATCH --job-name=copy_reads_parallel
#SBATCH --partition=batch                                # Partition (queue)
#SBATCH --ntasks=1                                       # Single task
#SBATCH --cpus-per-task=16                               # CPUs per task
#SBATCH --mem=40gb                                       # Memory
#SBATCH --time=05-00:00:00                               # Time limit (HH:MM:SS)
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out    # STDOUT log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err     # STDERR log
#SBATCH --mail-type=END,FAIL                             # Mail events
#SBATCH --mail-user=ma95362@uga.edu                      # Your email

module load parallel   # parallel *is* usually a module, keep this line

#SRC=/work/fdqlab/Ethiopia_dataset_763
SRC=/scratch/ma95362/eth_national_reads
#DEST=/scratch/ma95362/eth_national_reads
DEST=/scratch/ma95362/eth_national_analysis/all_fastq_reads
mkdir -p $DEST

cd $SRC || { echo "Failed to cd into $SRC"; exit 1; }

#find $SRC -type f | parallel -j 8 rsync -avh --partial {} $DEST/
find "$SRC" -type f -name "*.fastq.gz" | parallel -j 8 rsync -avh --partial {} "$DEST/"

