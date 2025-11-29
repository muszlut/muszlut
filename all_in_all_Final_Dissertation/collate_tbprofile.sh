#!/bin/bash
#SBATCH --job-name=collate_tbprofiler      
#SBATCH --partition=batch            
#SBATCH --ntasks=1                   
#SBATCH --cpus-per-task=8               
#SBATCH --mem=40gb                     
#SBATCH --time=06:00:00         
#SBATCH --output=/scratch/ma95362/all_in_all_reads/bactopia_prepare/my_tbprofiler_results/logs/collate/log.%j.out          
#SBATCH --error=/scratch/ma95362/all_in_all_reads/bactopia_prepare/my_tbprofiler_results/logs/collate/log.%j.err             
#SBATCH --mail-type=END,FAIL         
#SBATCH --mail-user=ma95362@uga.edu  
#Set output directory variable
OUTDIR="/scratch/ma95362/all_in_all_reads/bactopia_prepare/my_tbprofiler_results"
FOFN='/scratch/ma95362/all_in_all_reads/bactopia_prepare/my_tbprofiler_results'

#Load TB-Profiler module     

module load TB-Profiler/6.6.5
#Tell the program to make  the outdir folder
if [ ! -d $OUTDIR ] 
    then 
        mkdir -p $OUTDIR
fi

#Load modules
#module load Miniconda3/23.5.2-0

#activate conda envirnoment (should be version 6.3.0)

#source activate tb-profiler-env

#move to working directory:
cd $OUTDIR

#samples just needs to be a list of sample names. No path is required.
tb-profiler collate --dir $OUTDIR --samples samples_only.txt  --itol
