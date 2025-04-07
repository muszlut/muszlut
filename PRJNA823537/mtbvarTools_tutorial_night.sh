#!/bin/bash
#SBATCH --job-name=mtbvarTools_first_tutorial_job              # Job name
#SBATCH --partition=batch                                      # Partition (queue) name
#SBATCH --ntasks=1                                             # Run on a single CPU
#SBATCH --cpus-per-task=8                                      # Number of cores per task
#SBATCH --mem=40gb                                             # Job memory request
#SBATCH --time=07-00:00:00                                     # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err            # Standard error log

#SBATCH --mail-type=END,FAIL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                          # Where to send mail

#Set output directory variable
OUTDIR="/scratch/ma95362/PRJNA823537_ET125/mtbVARTolls_trial"

#Tell the program to make  the outdir folder
if [ ! -d $OUTDIR ] 
    then 
        mkdir -p $OUTDIR
fi

# Load required modules
module load fastp bwa samtools bcftools iqtree
source ~/mtbvartools_env/bin/activate  # Activate virtual environment

# Create output directory
cd $OUTDIR
# Step 1: Preprocessing FASTQ files
fastp -i /scratch/ma95362/PRJNA823537_ET125/ena-multiple-samples/fastqs/reads_R1.fastq \
      -I /scratch/ma95362/PRJNA823537_ET125/ena-multiple-samples/fastqs/reads_R2.fastq \
      -o /scratch/ma95362/PRJNA823537_ET125/mtbVARTolls_trial/trimmed_R1.fastq \
      -O /scratch/ma95362/PRJNA823537_ET125/mtbVARTolls_trial/trimmed_R2.fastq \
      -q 20 -l 50

# Step 2: Mapping reads to reference genome
bwa index reference.fasta
bwa mem reference.fasta /scratch/ma95362/PRJNA823537_ET125/mtbVARTolls_trial/trimmed_R1.fastq \
                        /scratch/ma95362/PRJNA823537_ET125/mtbVARTolls_trial/trimmed_R2.fastq \
                        > /scratch/ma95362/PRJNA823537_ET125/mtbVARTolls_trial/aligned_reads.sam
samtools view -Sb /scratch/ma95362/PRJNA823537_ET125/mtbVARTolls_trial/aligned_reads.sam | \
samtools sort -o /scratch/ma95362/PRJNA823537_ET125/mtbVARTolls_trial/aligned_reads.bam
samtools index /scratch/ma95362/PRJNA823537_ET125/mtbVARTolls_trial/aligned_reads.bam

# Step 3: Variant calling
bcftools mpileup -Ou -f reference.fasta /scratch/ma95362/PRJNA823537_ET125/mtbVARTolls_trial/aligned_reads.bam | \
bcftools call -mv -Oz -o /scratch/ma95362/PRJNA823537_ET125/mtbVARTolls_trial/variants.vcf.gz
bcftools index /scratch/ma95362/PRJNA823537_ET125/mtbVARTolls_trial/variants.vcf.gz

# Step 4: Phylogenetic analysis
bcftools consensus -f reference.fasta /scratch/ma95362/PRJNA823537_ET125/mtbVARTolls_trial/variants.vcf.gz > \
/scratch/ma95362/PRJNA823537_ET125/mtbVARTolls_trial/consensus.fasta
iqtree -s /scratch/ma95362/PRJNA823537_ET125/mtbVARTolls_trial/consensus.fasta -m GTR+G -bb 1000

echo "Pipeline completed successfully!"