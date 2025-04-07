#!/bin/bash
#SBATCH --job-name=mtbvarTools_batch_job              # Job name
#SBATCH --partition=batch                             # Partition (queue) name
#SBATCH --ntasks=1                                    # Run on a single CPU
#SBATCH --cpus-per-task=8                             # Number of cores per task
#SBATCH --mem=40gb                                    # Job memory request
#SBATCH --time=07-00:00:00                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out  # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err   # Standard error log

#SBATCH --mail-type=END,FAIL                          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu                   # Where to send mail

# Reference genome path
REFERENCE="/scratch/ma95362/reference/ncbi_dataset/data/GCF_000195955.2/genomic.fna"

# Input directory containing FASTQ files
INPUT_DIR="/scratch/ma95362/PRJNA823537_ET125/ena-multiple-samples/fastqs"

# Output directory for results
OUTDIR="/scratch/ma95362/PRJNA823537_ET125/mtbVARTolls_trial"

# Create output directory if it doesn't exist
if [ ! -d $OUTDIR ]; then
    mkdir -p $OUTDIR
fi

# Load required modules
module load fastp/0.23.2-GCC-11.3.0
module load BWA/0.7.17-GCCcore-12.3.0
module load SAMtools/1.18-GCC-12.3.0
module load BCFtools/1.18-GCC-12.3.0
source ~/mtbvartools_env/bin/activate  # Activate virtual environment

# Build BWA index (if not already built)
if [ ! -f ${REFERENCE}.bwt ]; then
    echo "Building BWA index for reference genome..."
    bwa index $REFERENCE
fi

# Change to output directory
cd $OUTDIR

# Loop over all FASTQ pairs
for SAMPLE in $(ls $INPUT_DIR/*_R1.fastq.gz | sed 's/_R1.fastq.gz//g'); do
    SAMPLE_NAME=$(basename $SAMPLE)
    echo "Processing sample: $SAMPLE_NAME"

    # Step 1: Preprocessing FASTQ files
    fastp -i ${SAMPLE}_R1.fastq.gz \
          -I ${SAMPLE}_R2.fastq.gz \
          -o $OUTDIR/${SAMPLE_NAME}_trimmed_R1.fastq.gz \
          -O $OUTDIR/${SAMPLE_NAME}_trimmed_R2.fastq.gz \
          -q 20 -l 50
    if [ $? -ne 0 ]; then
        echo "Error in fastp for $SAMPLE_NAME. Skipping..."
        continue
    fi

    # Step 2: Mapping reads to reference genome
    bwa mem $REFERENCE $OUTDIR/${SAMPLE_NAME}_trimmed_R1.fastq.gz $OUTDIR/${SAMPLE_NAME}_trimmed_R2.fastq.gz \
                        > $OUTDIR/${SAMPLE_NAME}_aligned_reads.sam
    if [ $? -ne 0 ]; then
        echo "Error in BWA MEM for $SAMPLE_NAME. Skipping..."
        continue
    fi

    samtools view -Sb $OUTDIR/${SAMPLE_NAME}_aligned_reads.sam | \
    samtools sort -o $OUTDIR/${SAMPLE_NAME}_aligned_reads.bam
    if [ $? -ne 0 ]; then
        echo "Error in SAMtools for $SAMPLE_NAME. Skipping..."
        continue
    fi

    samtools index $OUTDIR/${SAMPLE_NAME}_aligned_reads.bam
    if [ $? -ne 0 ]; then
        echo "Error in SAMtools indexing for $SAMPLE_NAME. Skipping..."
        continue
    fi

    # Step 3: Variant calling
    bcftools mpileup -Ou -f $REFERENCE $OUTDIR/${SAMPLE_NAME}_aligned_reads.bam | \
    bcftools call -mv -Oz -o $OUTDIR/${SAMPLE_NAME}_variants.vcf.gz
    if [ $? -ne 0 ]; then
        echo "Error in BCFtools variant calling for $SAMPLE_NAME. Skipping..."
        continue
    fi

    bcftools index $OUTDIR/${SAMPLE_NAME}_variants.vcf.gz
    if [ $? -ne 0 ]; then
        echo "Error in BCFtools indexing for $SAMPLE_NAME. Skipping..."
        continue
    fi

    echo "Finished processing sample: $SAMPLE_NAME"
done

echo "Batch processing completed for all samples!"