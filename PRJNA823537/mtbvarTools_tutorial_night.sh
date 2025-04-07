#!/bin/bash
#SBATCH --job-name=mtbvarTools_first_tutorial_job              # Job name
#SBATCH --partition=batch                                      # Partition (queue) name
#SBATCH --ntasks=1                                             # Run on a single CPU
#SBATCH --cpus-per-task=8                                      # Number of cores per task
#SBATCH --mem=40gb                                             # Job memory request
#SBATCH --time=07-00:00:00                                     # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err            # Standard error log

#SBATCH --mail-type=END,FAIL                                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                            # Where to send mail

# Set output directory variable
OUTDIR="/scratch/ma95362/PRJNA823537_ET125/mtbVARTolls_trial"

# Create the output directory if it doesn't exist
if [ ! -d $OUTDIR ]; then
    mkdir -p $OUTDIR
fi

# Load required modules
module load fastp/0.23.2-GCC-11.3.0
module load BWA/0.7.17-GCCcore-12.3.0
module load SAMtools/1.18-GCC-12.3.0
module load BCFtools/1.18-GCC-12.3.0
source ~/mtbvartools_env/bin/activate  # Activate virtual environment

# Change to the output directory
cd $OUTDIR

# Step 1: Preprocessing FASTQ files
fastp -i /scratch/ma95362/PRJNA823537_ET125/ena-multiple-samples/fastqs/reads_R1.fastq \
      -I /scratch/ma95362/PRJNA823537_ET125/ena-multiple-samples/fastqs/reads_R2.fastq \
      -o $OUTDIR/trimmed_R1.fastq \
      -O $OUTDIR/trimmed_R2.fastq \
      -q 20 -l 50

# Step 2: Mapping reads to reference genome
bwa index reference.fasta
bwa mem reference.fasta $OUTDIR/trimmed_R1.fastq $OUTDIR/trimmed_R2.fastq \
                        > $OUTDIR/aligned_reads.sam
samtools view -Sb $OUTDIR/aligned_reads.sam | \
samtools sort -o $OUTDIR/aligned_reads.bam
samtools index $OUTDIR/aligned_reads.bam

# Step 3: Variant calling
bcftools mpileup -Ou -f reference.fasta $OUTDIR/aligned_reads.bam | \
bcftools call -mv -Oz -o $OUTDIR/variants.vcf.gz
bcftools index $OUTDIR/variants.vcf.gz

echo "Pipeline completed successfully!"