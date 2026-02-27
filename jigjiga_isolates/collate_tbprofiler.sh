#!/bin/bash
#SBATCH --job-name=collate_tbprofiler       # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=8               # Number of cores per task
#SBATCH --mem=40gb                     # Job memory request
#SBATCH --time=06:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out           # Standard output log
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err             # Standard error log

#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ma95362@uga.edu  # Where to send mail	

#Set output directory variable
OUTDIR="/scratch/ma95362/publication/tbprofiler_results"
#FOFN='/scratch/ma95362/publication/tbprofiler_results'

# Load Conda/Miniforge module
module load Miniforge3

# Initialize Conda for bash
eval "$(conda shell.bash hook)"

# Activate TB-Profiler environment
conda activate tbprofiler_env

#Tell the program to make  the outdir folder
if [ ! -d $OUTDIR ] 
    then 
        mkdir -p $OUTDIR
fi
#move to working directory:
cd $OUTDIR
#to generate FOFN 
#cd /scratch/ma95362/publication/tbprofiler_results
#ls -d */ | sed 's#/##' > samples.txt

#samples just needs to be a list of sample names. No path is required.
#if FOFN is there
#tb-profiler collate --samples samples.txt --dir $OUTDIR --itol
#otherwise
tb-profiler collate \
    --dir /scratch/ma95362/publication/tbprofiler_results/results \
    --itol

#To isoltae the binary output of the spoligotype, you can use the following command in the /results directory directly on the terminarl or here :
#echo -e "Sample\tSpoligotype_binary" > spoligotype_binary_matrix.txt
#
#for f in *.results.json
#do
#    sample=$(basename $f .results.json)
#    binary=$(grep -A1 '"spoligotype"' $f | grep '"binary"' | head -1 | sed 's/.*"binary": "//;s/".*//')
#    echo -e "${sample}\t${binary}" >> spoligotype_binary_matrix.txt
#done

#AND, I used this to call all the spoligotypes in one file, but it is not formatted as a matrix. It is just a list of sample names and their spoligotypes.

#echo -e "SampleID\tLineageName\tLineageNumber\tSublineage\tBinary\tOctal\tFamily\tSIT" > tbprofiler_summary.tsv && \
#for f in *.txt; do \
#  id=$(grep -m1 "^ID:" "$f" | awk '{print $2}'); \
#  lineage_name=$(grep -m1 "^Lineage Fraction" -A1 "$f" | tail -n1 | awk '{print $3}'); \
#  lineage_number=$(grep -m1 "^Lineage Fraction" -A1 "$f" | tail -n1 | awk '{print $1}'); \
#  sublineage=$(grep -m1 "^Lineage Fraction" -A1 "$f" | tail -n1 | awk '{print $2}'); \
#  binary=$(grep -m1 "^Binary:" "$f" | awk '{print $2}'); \
#  octal=$(grep -m1 "^Octal:" "$f" | awk '{print $2}'); \
#  family=$(grep -m1 "^Family:" "$f" | awk '{print $2}'); \
#  sit=$(grep -m1 "^SIT:" "$f" | awk '{print $2}'); \
#  echo -e "${id}\t${lineage_name}\t${lineage_number}\t${sublineage}\t${binary}\t${octal}\t${family}\t${sit}"; \
#done >> tbprofiler_summary.tsv
