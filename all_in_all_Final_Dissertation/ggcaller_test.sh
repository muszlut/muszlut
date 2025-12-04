#!/bin/bash
#SBATCH --job-name=ggcaller_CPS
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/ma95362/ggcaller_module/logs/log.%j.out
#SBATCH --error=/scratch/ma95362/ggcaller_module/logs/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

# ==============================================================================
# ሞጁሎችን መጫን (Module Loading)
# ==============================================================================
# ggCaller እና ሌሎች የሚያስፈልጉ ሶፍትዌሮችን (ለምሳሌ Python, Panaroo) መጫን።
# ክላስተርዎ ሞጁል ሲስተም የሚጠቀም ከሆነ፣ እንደ ክላስተሩ ህግ እነዚህን ማስተካከል ይኖርብዎታል።
module load ggCaller/1.4.1

# ==============================================================================
# የሥራ ቦታ (Working Directory)
# ==============================================================================
# ወደ ፋይሎችዎ ቦታ መሄድ
cd /scratch/ma95362/all_in_all_reads/test_reads

# ------------------------------------------------------------------------------
# ደረጃ 1: ለ ggCaller የሚያስፈልገውን የፋይል ዝርዝር ማዘጋጀት
# ------------------------------------------------------------------------------

echo "1. የፋይል ዝርዝር እየተዘጋጀ ነው (input.txt)..."

# $PWD ማለት አሁን ያለንበት ፎልደር አድራሻ ማለት ነው።
# .fastq.gz የሚጨርሱትን ፋይሎች በሙሉ ይዘረዝራል።
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
# ደረጃ 2: ggCallerን ማስኬድ
# ------------------------------------------------------------------------------

echo "2. ggCaller በ --reads mode እየተጀመረ ነው..."

# $SLURM_CPUS_PER_TASK ከላይ በ #SBATCH የተሰጠውን የኮር ብዛት (8) ይጠቀማል።
ggcaller --reads input.txt \
         --threads 8 \
         --out Reads_Pangenome_Output \
         --balrog-db /scratch/ma95362/ggcaller_db/ggCallerdb \
         --annotation sensitive # ጂኖችን በጥልቀት ስም ለመስጠት ይህንን አማራጭ እንጨምራለን

# ------------------------------------------------------------------------------
# መጨረሻ
# ------------------------------------------------------------------------------

if [ $? -eq 0 ]; then
    echo "3. ggCaller ሥራውን በተሳካ ሁኔታ አጠናቋል።"
    echo "ውጤቶችዎ በ Reads_Pangenome_Output ፎልደር ውስጥ ይገኛሉ።"
else
    echo "4. !!! ስህተት: ggCaller ሥራውን ሳያጠናቅቅ ቆሟል። እባክዎ slurm-*.err ፋይሉን ይመልከቱ።"
fi