#!/bin/bash
#SBATCH --job-name=CollateTBprofiler
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40gb
#SBATCH --time=07-00:00:00
#SBATCH --output=/scratch/ma95362/scratch/log.%j.out
#SBATCH --error=/scratch/ma95362/scratch/log.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ma95362@uga.edu

set -euo pipefail

# Load Micromamba
module load Micromamba/2.3.0
eval "$(micromamba shell hook --shell bash)"
micromamba activate tbprofiler

# Set directories
OUTDIR="/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples/all_fastqs/TBprofiler_results_conda"
FOFN="/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples/all_fastqs/tbprofiler.fofn"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "Running TB-Profiler collation..."
tb-profiler collate \
  --file_list "$FOFN" \
  --prefix selected_tbprofiler_results \
  --itol

# Generate lineage and sublineage frequency summary
echo "Summarizing lineage distribution..."
cut -f2,3 selected_tbprofiler_results.metadata.txt | tail -n +2 | \
  awk -F'\t' '{count[$1]++} END {for (i in count) print i, count[i]}' | sort -k2 -nr > lineage_counts.txt

cut -f3 selected_tbprofiler_results.metadata.txt | tail -n +2 | \
  awk '{count[$1]++} END {for (i in count) print i, count[i]}' | sort -k2 -nr > sublineage_counts.txt

# Calculate percentages
TOTAL=$(($(wc -l < selected_tbprofiler_results.metadata.txt)-1))
awk -v total=$TOTAL '{printf "%s\t%.2f%%\n", $0, ($2/total)*100}' lineage_counts.txt > lineage_percentage.txt
awk -v total=$TOTAL '{printf "%s\t%.2f%%\n", $0, ($2/total)*100}' sublineage_counts.txt > sublineage_percentage.txt

# Combine spoligotyping results (if TB-Profiler version ≥ 4.0)
if grep -q "spoligotype" selected_tbprofiler_results.metadata.txt; then
  echo "Extracting spoligotyping info..."
  cut -f1,15,16,17 selected_tbprofiler_results.metadata.txt > spoligotyping_summary.txt
  # (adjust columns depending on actual spoligotype, SIT, and octal column positions)
else
  echo "⚠️ Spoligotyping data not found in metadata — ensure TB-Profiler was run with --spoligotype flag."
fi

micromamba deactivate

echo "Done! Output files:"
echo "- selected_tbprofiler_results.metadata.txt"
echo "- lineage_percentage.txt"
echo "- sublineage_percentage.txt"
echo "- spoligotyping_summary.txt (if available)"
