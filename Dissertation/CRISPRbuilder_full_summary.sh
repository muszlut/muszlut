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

# Header
echo -e "sample_id\tspacer_count\tgaps\treads_with_2_spacers\tspacer_order\tbinary_43\tdeleted_spacers" \
> CRISPRbuilder_full_summary.tsv

# Loop through isolate directories
for d in */ ; do
    sample=$(basename "$d")

    pattern_file="${d}/${sample}_crispr_patterns.blast"
    reads_file="${d}/${sample}.reads_with_2_spacers"
    gaps_file="${d}/${sample}.not_consecutive"

    spacer_order="NA"
    spacer_count="0"
    reads="0"
    gaps="0"
    binary_43=""
    deleted=""

    # Reads with ≥2 spacers (confidence)
    [[ -f "$reads_file" ]] && reads=$(wc -l < "$reads_file")

    # CRISPR gaps
    [[ -f "$gaps_file" ]] && gaps=$(wc -l < "$gaps_file")

    # Extract ordered unique spacers
    if [[ -f "$pattern_file" ]]; then
        spacer_order=$(awk '{print $2}' "$pattern_file" | uniq | tr '\n' ',')
        spacer_count=$(echo "$spacer_order" | awk -F',' '{print NF-1}')

        # Build binary 43-spacer format
        for i in $(seq 1 43); do
            if echo "$spacer_order" | grep -qw "$i"; then
                binary_43="${binary_43}1"
            else
                binary_43="${binary_43}0"
                deleted="${deleted}${i},"
            fi
        done
    fi

    # Clean trailing commas
    spacer_order=${spacer_order%,}
    deleted=${deleted%,}

    echo -e "${sample}\t${spacer_count}\t${gaps}\t${reads}\t${spacer_order}\t${binary_43}\t${deleted}" \
    >> CRISPRbuilder_full_summary.tsv
done
