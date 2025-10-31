from panaroo.plotting import rarefaction_plot

# Define input and output paths
input_dir = "/scratch/ma95362/eth_national_analysis/all_fastq_reads/ETH_paired_end_samples/bactopia-runs/pangenome_of_1368/panaroo"
output_file = input_dir + "/Plot_output/panaroo_rarefaction_curve.png"

# Generate the rarefaction plot
rarefaction_plot(input_dir=input_dir, output_file=output_file)