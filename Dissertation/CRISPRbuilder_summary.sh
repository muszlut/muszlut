#!/bin/bash
#SBATCH --job-name=prep_shuffled_fasta
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/ma95362/CRISPRbuilder-TB/logs/logs.%A_%a.out
#SBATCH --error=/scratch/ma95362/CRISPRbuilder-TB/logs/logs.%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu
cd /scratch/ma95362/CRISPRbuilder-TB/sequences

echo -e "sample_id\tspacer_count\tgaps\treads_with_2_spacers\tcrispr_pattern_file" > CRISPRbuilder_summary.tsv

for d in */ ; do
    sample=$(basename "$d")

    pattern_file="${d}/${sample}_crispr_patterns.blast"
    reads_file="${d}/${sample}.reads_with_2_spacers"
    gaps_file="${d}/${sample}.not_consecutive"

    spacer_count="NA"
    reads="0"
    gaps="0"

    if [[ -f "$pattern_file" ]]; then
        spacer_count=$(grep -c "Spacer" "$pattern_file")
    fi

    if [[ -f "$reads_file" ]]; then
        reads=$(wc -l < "$reads_file")
    fi

    if [[ -f "$gaps_file" ]]; then
        gaps=$(wc -l < "$gaps_file")
    fi

    echo -e "${sample}\t${spacer_count}\t${gaps}\t${reads}\t${pattern_file}" >> CRISPRbuilder_summary.tsv
done
