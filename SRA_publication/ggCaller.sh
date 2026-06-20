#!/bin/bash
#SBATCH --job-name=ggCaller
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120gb
#SBATCH --time=07-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# ==============================================================================
# የሥራ ቦታ (Working Directory)
# ==============================================================================
# አሁን ባሉበት ፎልደር መሰረት ተስተካክሏል
cd /scratch/ma95362/257_assembled_files

# ------------------------------------------------------------------------------
# ደረጃ 1: ለ ggCaller የሚያስፈልገውን የፋይል ዝርዝር ማዘጋጀት
# ------------------------------------------------------------------------------

#echo "1. የተገጣጠሙ ፋይሎች (Assembly FASTA/fna) ዝርዝር እየተዘጋጀ ነው (input.txt)..."

# ማስተካከያ: .fasta ወደ .fna.gz ተቀይሯል. ggCaller የ .gz (Gzip) ቅርጸትን ማንበብ ይችላል።
#ls -d -1 $PWD/*.fna > input.txt 

#if [ -f "input.txt" ]; then
#    echo "ዝርዝሩ በ input.txt ላይ ተቀምጧል:"
#    wc -l input.txt
#else
#    echo "!!! ስህተት: input.txt አልተፈጠረም። የፋይል ቅጥያዎ (.fna) ትክክል መሆኑን ያረጋግጡ።"
#    exit 1
#fi

# ------------------------------------------------------------------------------
# ደረጃ 2: ggCallerን ማስኬድ (--refs ተጠቅመን)
# ------------------------------------------------------------------------------

#echo "2. ggCaller በ --refs mode (Assembly FNA) እየተጀመረ ነው..."

# --refs የሚለው አማራጭ Assembly ፋይሎችን ማንበብ እንዲችል ያደርገዋል።
#singularity exec /apps/singularity-images/ggcaller_v1.5.0.sif \
#ggcaller \
#  --refs input.txt \
#  --gene-finding-only \
#  --save \
#  --threads 32 \
#  --balrog-db /scratch/ma95362/ggcaller_db \
#  --out MTB_ggcaller


# ------------------------------------------------------------------------------
# መጨረሻ
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ደረጃ 3: ggCallerን ማስኬድ (--refs ተጠቅመን)
# ------------------------------------------------------------------------------

#echo "3. ggCaller በ core genome alignment እየተጀመረ ነው..."

singularity exec /apps/singularity-images/ggcaller_v1.5.0.sif \
ggcaller \
  --refs input.txt \
  --save \
  --threads 32 \
  --balrog-db /scratch/ma95362/ggcaller_db \
  --alignment core \
  --aligner def \
  --core-threshold 0.99 \
  --out Alignment_results