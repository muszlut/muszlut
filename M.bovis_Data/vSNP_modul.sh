#!/bin/bash
#SBATCH --job-name=vSNP_bovis
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Set output directory and create if it doesn't exist
OUTDIR="/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples/all_fastqs"
REF_DIR="/home/ma95362/vsnp3_test_dataset/vsnp_dependencies/Mycobacterium_AF2122"
Step1="/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples/all_fastqs/vSNP_output"/step1_output
Step2="/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples/all_fastqs/vSNP_output"/step2_output
REF_FASTA=${REF_DIR}/NC_002945v4.fasta      
mkdir -p $OUTDIR $Step1 $Step2

module load vsnp3/3.26

cd $OUTDIR


# Run vSNP step 1 for each sample
for R1_file in $OUTDIR/*_R1.fastq.gz; do
    R2_file="${R1_file/_R1/_R2}"
    vsnp3_step1.py -r1 $R1_file -r2 $R2_file -f $REF_FASTA -o "$Step1" -spoligo
done

#Find all vcf files in the step1 output directory and move them in step2 directory
find $Step1 -type f -name "*_zc.vcf" -exec cp {} $Step2 \;
#move to step2 folder
cd $Step2
# Run vSNP step 2 to combine SNPs and build tree        
vsnp3_step2.py -a -t $REF_FASTA -remove_by_name $OUTDIR/bactopia_exclude.tsv