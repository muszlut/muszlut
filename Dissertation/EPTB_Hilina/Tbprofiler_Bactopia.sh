#!/bin/bash
#SBATCH --job-name=tbprofiler_by_module       
#SBATCH --partition=batch            
#SBATCH --ntasks=1                   
#SBATCH --cpus-per-task=16              
#SBATCH --mem=40gb                    
#SBATCH --time=06:00:00              
#SBATCH --output=/scratch/ma95362/EPTB_Hilina/new_project_logs/curl_download_all/Bactopia_prepare/log.%j.out          
#SBATCH --error=/scratch/ma95362/EPTB_Hilina/new_project_logs/curl_download_all/Bactopia_prepare/log.%j.err             

#SBATCH --mail-type=END,FAIL          
#SBATCH --mail-user=ma95362@uga.edu  
#Set output directory variable
OUTDIR="/scratch/ma95362/EPTB_Hilina/new_project_logs/curl_download_all/Bactopia_prepare/"
#Tell the program to make  the outdir folder
if [ ! -d $OUTDIR ] 
    then 
        mkdir -p $OUTDIR
fi
module load Bactopia/3.2.0-conda
#module load EDirect/20.5.20231006-GCCcore-12.3.0
#module load BBMap/39.19-GCC-13.3.0
cd $OUTDIR
bactopia \
    --wf tbprofiler \
    --bactopia $OUTDIR \
    --max_cpus 16
