#!/bin/bash
#SBATCH --job-name=ggcaller_CPS           # የስራው ስም
#SBATCH --partition=batch                 # የሚላክበት ክፍል (Partition)
#SBATCH --ntasks=1                        # ዋናው ፕሮሰስ ብዛት
#SBATCH --cpus-per-task=8               # የሚያስፈልግ የኮር ብዛት
#SBATCH --mem=32G                       # የሚያስፈልገው የማስታወሻ መጠን
#SBATCH --time=24:00:00                   # ከፍተኛው የሥራ ጊዜ (24 ሰዓት)
#SBATCH --output=/scratch/ma95362/ggcaller_module/logs/log.%j.out
#SBATCH --error=/scratch/ma95362/ggcaller_module/logs/log.%j.err
#SBATCH --mail-type=END,FAIL              # ስራው ሲያልቅ ወይም ሲበላሽ ኢሜይል ይላካል
#SBATCH --mail-user=ma95362@uga.edu       # ኢሜይል የሚላክለት አድራሻ

# ==============================================================================
# ሞጁሎችን መጫን (Module Loading)
# ==============================================================================
# እርስዎ የመረጡት የ ggCaller ሞጁል ስሪት
module load ggCaller/1.4.1
# ggCaller/Panaroo የሚጠቀምባቸውን ሌሎች ሶፍትዌሮች ሊያካትት ይችላል
# (እርስዎ anaconda3ን ስላስወገዱ፣ ggCaller/1.4.1 ሁሉንም ጥገኞች እንደሚጭን ተገምቷል)

# ==============================================================================
# የሥራ ቦታ (Working Directory)
# ==============================================================================
# ወደ ፋይሎችዎ ቦታ መሄድ
cd /scratch/ma95362/all_in_all_reads/test_reads

# ------------------------------------------------------------------------------
# ደረጃ 1: ለ ggCaller የሚያስፈልገውን የፋይል ዝርዝር ማዘጋጀት
# ------------------------------------------------------------------------------

echo "1. የፋይል ዝርዝር እየተዘጋጀ ነው (input.txt)..."

# አሁን ያለንበትን ፎልደር (Working Directory) ተጠቅሞ የፋይል ዝርዝር ይፈጥራል።
ls -d -1 $PWD/*.fastq.gz > input.txt

if [ -f "input.txt" ]; then
    echo "ዝርዝሩ በ input.txt ላይ ተቀምጧል:"
    # የተዘረዘሩትን ፋይሎች ብዛት ማሳየት (ለማጣራት)
    wc -l input.txt
else
    echo "!!! ስህተት: input.txt አልተፈጠረም። የፋይል ቅጥያዎ ትክክል መሆኑን ያረጋግጡ።"
    exit 1
fi

# ------------------------------------------------------------------------------
# ደረጃ 2: ggCallerን ማስኬድ (ችግሩን የሚፈታው --clean-mode sensitive ተጨምሯል)
# ------------------------------------------------------------------------------

echo "2. ggCaller በ --reads mode እና Sensitive Clean Mode እየተጀመረ ነው..."

# --clean-mode sensitive የተጨመረው "No genes detected in graph" የሚለውን ችግር ለመፍታት ነው።
ggcaller --reads input.txt \
         --threads 8 \
         --out Reads_Pangenome_Output \
         --balrog-db /scratch/ma95362/ggcaller_db/ggCallerdb \
         --annotation sensitive \
         --clean-mode sensitive 

# ------------------------------------------------------------------------------
# መጨረሻ
# ------------------------------------------------------------------------------

if [ $? -eq 0 ]; then
    echo "3. ggCaller ሥራውን በተሳካ ሁኔታ አጠናቋል።"
    echo "ውጤቶችዎ በ Reads_Pangenome_Output ፎልደር ውስጥ ይገኛሉ።"
else
    echo "4. !!! ስህተት: ggCaller ሥራውን ሳያጠናቅቅ ቆሟል። እባክዎ slurm-*.err ፋይሉን ይመልከቱ።"
fi