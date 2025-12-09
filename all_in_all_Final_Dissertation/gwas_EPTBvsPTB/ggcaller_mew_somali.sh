#!/bin/bash
#SBATCH --job-name=ggcaller_newsomali_only          
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=895G
#SBATCH --time=7-00:00:00
#SBATCH --output=/scratch/ma95362/ggcaller_module/logs/somali/log.%j.out
#SBATCH --error=/scratch/ma95362/ggcaller_module/logs/somali/log.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL             # ስራው ሲያልቅ ወይም ሲበላሽ ኢሜይል ይላካል
#SBATCH --mail-user=ma95362@uga.edu       # ኢሜይል የሚላክለት አድራሻ
# ==============================================================================
# ሞጁሎችን መጫን (Module Loading)
# ==============================================================================
module load ggCaller/1.4.1

# ==============================================================================
# የሥራ ቦታ (Working Directory)
# ==============================================================================
# ggCaller የሚሰራበት ፎልደር
cd /scratch/ma95362/all_in_all_reads/ggcaller_pangenome_fna

# ------------------------------------------------------------------------------
# ደረጃ 1: ለ ggCaller የሚያስፈልገውን የፋይል ዝርዝር ማዘጋጀት
# ------------------------------------------------------------------------------

echo "1. የተገጣጠሙ ፋይሎች (Assembly FASTA/FNA) ዝርዝር እየተዘጋጀ ነው (input.txt)..."

# .fna.gz የሚጨርሱትን የ Assembly ፋይሎች ይዘረዝራል.
ls -d -1 $PWD/*.fna.gz > input.txt 

if [ [ -f "input.txt" ] ]; then
 echo "ዝርዝሩ በ input.txt ላይ ተቀምጧል:"
 wc -l input.txt
else
echo "!!! ስህተት: input.txt አልተፈጠረም። የፋይል ቅጥያዎ (.fna.gz) ትክክል መሆኑን ያረጋግጡ።"
exit 1
fi

# ------------------------------------------------------------------------------
# ደረጃ 2: ggCallerን ማስኬድ (በተጠቃሚው ጥብቅ መለኪያዎች)
# ------------------------------------------------------------------------------

echo "2. ggCaller በከፍተኛ ጥብቅነት (Strict Filtering) እየተጀመረ ነው..."

# የተጠቃሚው ጥያቄ መለኪያዎች በሙሉ ተካተዋል:
ggcaller --genomes input.txt \
        --threads 8 \
        --out Assembly_Pangenome_Output_strict \
        --min-gene-len 180 \
        --min-hmm-score 60 \
        --min-protein-length 70 \
        --cds-only \
        --keep-paralogs \
        --cluster                     
# ------------------------------------------------------------------------------
# መጨረሻ
# ------------------------------------------------------------------------------

if [ $? -eq 0 ]; then
 echo "3. ggCaller ሥራውን በተሳካ ሁኔታ አጠናቋል።"
  echo "ውጤቶችዎ በ Assembly_Pangenome_Output_strict ፎልደር ውስጥ ይገኛሉ።"
else
echo "4. !!! ስህተት: ggCaller ሥራውን ሳያጠናቅቅ ቆሟል። እባክዎ log_strict_v8.*.err ፋይሉን ይመልከቱ።"
fi
