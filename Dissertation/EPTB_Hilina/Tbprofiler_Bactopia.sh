#!/bin/bash
#SBATCH --job-name=bactopia_singularity       
#SBATCH --partition=batch            
#SBATCH --ntasks=1                   
#SBATCH --cpus-per-task=8              
#SBATCH --mem=40gb                    
#SBATCH --time=06:00:00              
#SBATCH --output=/scratch/ma95362/EPTB_Hilina/Bactopia_Run/split_reads/log.%j.out          
#SBATCH --error=/scratch/ma95362/EPTB_Hilina/Bactopia_Run/split_reads/log.%j.err             

#SBATCH --mail-type=END,FAIL          
#SBATCH --mail-user=ma95362@uga.edu  
#Set output directory variable
OUTDIR="/scratch/ma95362/EPTB_Hilina/Bactopia_Run/split_reads/"
#Tell the program to make  the outdir folder
if [ ! -d $OUTDIR ] 
    then 
        mkdir -p $OUTDIR
fi
module load Bactopia/3.2.0-conda
module load EDirect/20.5.20231006-GCCcore-12.3.0
cd $OUTDIR
bactopia \
    -profile singularity \
    --wf tbprofiler \
    --bactopia $OUTDIR
