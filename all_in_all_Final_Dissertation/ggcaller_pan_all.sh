#!/bin/bash
#SBATCH --job-name=ggcaller_all           # የስራው ስም
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=2000G
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/ma95362/ggcaller_module/logs/log.%j.out
#SBATCH --error=/scratch/ma95362/ggcaller_module/logs/log.%j.err
#SBATCH --mail-type=END,FAIL              # ስራው ሲያልቅ ወይም ሲበላሽ ኢሜይል ይላካል
#SBATCH --mail-user=ma95362@uga.edu       # ኢሜይል የሚላክለት አድራሻ

# ==============================================================================
# ሞጁሎችን መጫን (Module Loading)
# ==============================================================================
module load ggCaller/1.4.1

# ==============================================================================
# የሥራ ቦታ (Working Directory)
# ==============================================================================
# አሁን ባሉበት ፎልደር መሰረት ተስተካክሏል
cd /scratch/ma95362/all_in_all_reads/ggcaller_pangenome_fna

# ------------------------------------------------------------------------------
# ደረጃ 1: ለ ggCaller የሚያስፈልገውን የፋይል ዝርዝር ማዘጋጀት
# ------------------------------------------------------------------------------

echo "1. የተገጣጠሙ ፋይሎች (Assembly FASTA/FNA) ዝርዝር እየተዘጋጀ ነው (input.txt)..."

# ማስተካከያ: .fasta ወደ .fna.gz ተቀይሯል. ggCaller የ .gz (Gzip) ቅርጸትን ማንበብ ይችላል።
ls -d -1 $PWD/*.fna.gz > input.txt 

if [ -f "input.txt" ]; then
    echo "ዝርዝሩ በ input.txt ላይ ተቀምጧል:"
    wc -l input.txt
else
    echo "!!! ስህተት: input.txt አልተፈጠረም። የፋይል ቅጥያዎ (.fna.gz) ትክክል መሆኑን ያረጋግጡ።"
    exit 1
fi

# ------------------------------------------------------------------------------
# ደረጃ 2: ggCallerን ማስኬድ (--refs ተጠቅመን)
# ------------------------------------------------------------------------------

echo "2. ggCaller በ --refs mode (Assembly FNA) እየተጀመረ ነው..."

# --refs የሚለው አማራጭ Assembly ፋይሎችን ማንበብ እንዲችል ያደርገዋል።
ggcaller --refs input.txt \
         --threads 8 \
         --out Assembly_Pangenome_Output \
         --balrog-db /scratch/ma95362/ggcaller_db/ggCallerdb \
         --annotation sensitive 

# ------------------------------------------------------------------------------
# መጨረሻ
# ------------------------------------------------------------------------------

if [ $? -eq 0 ]; then
    echo "3. ggCaller ሥራውን በተሳካ ሁኔታ አጠናቋል።"
    echo "ውጤቶችዎ በ Assembly_Pangenome_Output ፎልደር ውስጥ ይገኛሉ።"
else
    echo "4. !!! ስህተት: ggCaller ሥራውን ሳያጠናቅቅ ቆሟል። እባክዎ slurm-*.err ፋይሉን ይመልከቱ።"
fi