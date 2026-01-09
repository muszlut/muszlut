#!/bin/bash
#SBATCH --job-name=prep_crisprbuilder_spoligo
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

# ============================================================
# CRISPRbuilder → Classical Spoligotype Converter
# Author: Musse
# ============================================================

# ---- INPUTS ----
CRISPR_DIR="/scratch/ma95362/CRISPRbuilder-TB/sequences"      # directory with sample folders
SIT_MAP="/scratch/ma95362/CRISPRbuilder-TB/data/SIT.xls"      # CRISPRbuilder spacer → classical spacer
OUT="/scratch/ma95362/CRISPRbuilder-TB/CRISPRbuilder_spoligotype_table.tsv"

# ---- HEADER ----
echo -e "sample_id\tcb_spacers\tclassical_spacers\tbinary_43\toctal\tdeleted_spacers" > $OUT

# ---- LOOP THROUGH SAMPLES ----
for sdir in "$CRISPR_DIR"/*/; do
    sample=$(basename "$sdir")
    blast_file="$sdir/${sample}_crispr_patterns.blast"
    [[ ! -f "$blast_file" ]] && { echo "⚠️  Blast file not found for $sample, skipping"; continue; }

    echo "Processing $sample ..."

    # --------------------------------------------------------
    # 1. Extract CRISPRbuilder esp spacers (last two columns)
    # --------------------------------------------------------
    awk -F',' '{print $(NF-1); print $NF}' "$blast_file" \
        | sort -n | uniq > "$sdir/${sample}.cb_spacers"

    cb_spacers=$(paste -sd, "$sdir/${sample}.cb_spacers")

    # --------------------------------------------------------
    # 2. Map CRISPRbuilder → classical spacers (1–43)
    # --------------------------------------------------------
    > "$sdir/${sample}.classical_spacers"

    while read cb; do
        awk -F'\t' -v cb="$cb" '$2==cb {print $1}' "$SIT_MAP" \
            >> "$sdir/${sample}.classical_spacers"
    done < "$sdir/${sample}.cb_spacers"

    sort -n "$sdir/${sample}.classical_spacers" | uniq > "$sdir/tmp"
    mv "$sdir/tmp" "$sdir/${sample}.classical_spacers"

    classical=$(paste -sd, "$sdir/${sample}.classical_spacers")

    # --------------------------------------------------------
    # 3. Build binary-43 + deletion list
    # --------------------------------------------------------
    bin=""
    del=""

    for i in $(seq 1 43); do
        if grep -qw "$i" "$sdir/${sample}.classical_spacers"; then
            bin="${bin}1"
        else
            bin="${bin}0"
            del="${del}${i},"
        fi
    done
    del=${del%,}  # remove trailing comma

    # --------------------------------------------------------
    # 4. Convert binary-43 → octal
    # --------------------------------------------------------
    octal=""
    for i in $(seq 0 14); do
        chunk=${bin:$((i*3)):3}
        [[ ${#chunk} -lt 3 ]] && chunk=$(printf "%-3s" "$chunk" | tr ' ' '0')
        octal="${octal}$(echo "ibase=2; $chunk" | bc)"
    done

    # --------------------------------------------------------
    # 5. Write output row
    # --------------------------------------------------------
    echo -e "$sample\t$cb_spacers\t$classical\t$bin\t$octal\t$del" >> $OUT

done

echo "DONE ✔  Output written to: $OUT"
