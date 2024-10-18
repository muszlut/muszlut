#Set output directory variable
OUTDIR="/scratch/ma95362/gbk"
SCRIPT="/home/ma95362/muszlut"


#Load modules 
module load Biopython/1.84-foss-2023b


#move to working diectory
cd $OUTDIR

Python $SCRIPT/convert_gbff_to_gbk.py $OUTDIR/genomic.gbff $OUTDIR/genomic.gbk
