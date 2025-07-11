import os
import yaml

# CONFIGURE THESE:
study_name = "Mbovis_Study"
output_dir = "/scratch/ma95362/ETH_bovis_Sequence/magma_output"
reads_dir = "/scratch/ma95362/ETH_bovis_Sequence/reads/64_SH_Sequence_data/raw/sub_raw"
yaml_out = "my_parameters_1.yml"

# find all R1 files
r1_files = sorted([f for f in os.listdir(reads_dir) if f.endswith("_R1.fastq.gz")])

samples = []

for r1 in r1_files:
    sample_id = r1.replace("_R1.fastq.gz", "")
    r2 = r1.replace("_R1.fastq.gz", "_R2.fastq.gz")
    r1_path = os.path.join(reads_dir, r1)
    r2_path = os.path.join(reads_dir, r2)

    # sanity check: does R2 exist?
    if not os.path.exists(r2_path):
        print(f"⚠️  Warning: R2 file missing for {sample_id}")
        continue

    samples.append({
        'sample_id': sample_id,
        'library': 1,
        'attempt': 1,
        'r1': r1_path,
        'r2': r2_path,
        'flowcell': 1,
        'lane': 1,
        'index_sequence': 1
    })

data = {
    'study_name': study_name,
    'output_dir': output_dir,
    'samples': samples
}

with open(yaml_out, 'w') as f:
    yaml.dump(data, f, sort_keys=False)

print(f"✅ YAML file written to: {yaml_out}")
