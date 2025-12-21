#!/bin/bash
#SBATCH --job-name=ncbi_upload      
#SBATCH --partition=batch           
#SBATCH --ntasks=1                  
#SBATCH --mem=4gb                   
#SBATCH --time=48:00:00             
#SBATCH --output=ncbi_upload_%j.out
#SBATCH --error=ncbi_upload_%j.err

# 1. Load the Miniforge module
ml Miniforge3/24.7.1-0

# 2. Activate your environment
source activate ~/ncbi_upload_env

# 3. Your Specific NCBI Credentials
NCBI_USER="subftp"
NCBI_PASS="CrasHeshyevVafivgud2"
REMOTE_DIR="uploads/musse.girma_aau.edu.et_vf8TknxA"

# 4. Give your subfolder a meaningful name (NO SPACES)
# You can change 'my_fastq_submission' to something like 'Girma_Project_2024'
SUBMISSION_FOLDER="Ethiopian_Somali_region_TB_fastq_submission"

# 5. Navigate to your data directory on Sapelo2
# Replace the path below with the actual folder where your FASTQs are
cd /scratch/ma95362/YOUR_DATA_FOLDER_NAME

echo "Starting upload at $(date)"

# 6. Perform the transfer using lftp
# We use 'mkdir -p' to create the subfolder and 'cd' to enter it
lftp -u ${NCBI_USER},${NCBI_PASS} ftp-private.ncbi.nlm.nih.gov <<EOF
set ftp:ssl-allow no
cd ${REMOTE_DIR}
mkdir -p ${SUBMISSION_FOLDER}
cd ${SUBMISSION_FOLDER}
mput *.fastq.gz
bye
EOF

echo "Upload completed at $(date)"