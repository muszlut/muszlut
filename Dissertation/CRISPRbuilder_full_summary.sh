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
#!/bin/bash

SIT_MAP="/scratch/ma95362/CRISPRbuilder-TB/data/SIT.xls"
cd /scratch/ma95362/CRISPRbuilder-TB/sequences
echo -e "sample_id\tcrisprbuilder_spacers\tcb_spacer_count\tbinary_43\tdeleted_spacers" \
> CRISPRbuilder_full_summary.tsv

for d in */ ; do
    sample=$(basename "$d")
    pattern_file="${d}/${sample}_crispr_patterns.blast"

    cb_spacers=""
    cb_count=0
    binary_43=""
    deleted=""

    if [[ -f "$pattern_file" ]]; then
        # Extract CRISPRbuilder spacer IDs (last two comma-separated fields)
        cb_spacers=$(awk -F',' '{print $(NF-1); print $NF}' "$pattern_file" \
                     | sort -n | uniq)

        cb_count=$(echo "$cb_spacers" | wc -l)

        # Build binary 43-spacer profile using SIT.xls
        for i in $(seq 1 43); do
            # Column 2 in SIT.xls = CRISPRbuilder spacer ID
            mapped=$(awk -F'\t' -v s="$i" '$1==s {print $2}' "$SIT_MAP")

            if echo "$cb_spacers" | grep -qw "$mapped"; then
                binary_43="${binary_43}1"
            else
                binary_43="${binary_43}0"
                deleted="${deleted}${i},"
            fi
        done
    fi

    deleted=${deleted%,}
    cb_spacers=$(echo "$cb_spacers" | tr '\n' ',' | sed 's/,$//')

    echo -e "${sample}\t${cb_spacers}\t${cb_count}\t${binary_43}\t${deleted}" \
    >> CRISPRbuilder_full_summary.tsv
done
