#!/bin/bash
#SBATCH --job-name=Tbprofiler_jason_collate_conversion
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40gb
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# Set variables
OUTDIR="/scratch/ma95362/PRJNA823537_ET125/Tbprofiler_module/results"
SCRIPT="/home/ma95362/muszlut/PRJNA823537/Tbprofiler_collate_from_new_module.py"

# Load modules
module load Biopython/1.84-foss-2023b

# Optional: ensure tqdm is available
pip install --user tqdm

# Move to output directory
cd $OUTDIR

# Run the script
python $SCRIPT --dir $OUTDIR
