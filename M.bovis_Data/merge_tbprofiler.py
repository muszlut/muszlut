import os
import pandas as pd

# Set the root directory where your sample folders are
root_dir = "/scratch/ma95362/ETH_M.bovis/m.bovis_Bactopia_Analysis/with_fixed_reads/M.bovis_paired_end_samples/all_fastqs/TBprofiler_results_conda"

# Create a list to hold dataframes
df_list = []

# Walk through all sample folders
for sample_folder in os.listdir(root_dir):
    sample_path = os.path.join(root_dir, sample_folder, "results")
    if os.path.exists(sample_path):
        # Find the TXT result file
        for file in os.listdir(sample_path):
            if file.endswith(".results.txt"):
                file_path = os.path.join(sample_path, file)
                # Read the TXT file into a dataframe
                df = pd.read_csv(file_path, sep="\t")
                df["Sample"] = sample_folder  # add sample column
                df_list.append(df)

# Concatenate all dataframes
merged_df = pd.concat(df_list, ignore_index=True)

# Save the merged table
output_file = os.path.join(root_dir, "TBProfiler_all_samples_merged.csv")
merged_df.to_csv(output_file, index=False)

print(f"Merged results saved to {output_file}")
