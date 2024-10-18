import sys
from Bio import SeqIO

# Get the input and output file paths from command-line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# Parse the .gbff file and write it as a .gbk file
with open(input_file, "r") as gbff_file, open(output_file, "w") as gbk_file:
    sequences = SeqIO.parse(gbff_file, "genbank")
    SeqIO.write(sequences, gbk_file, "genbank")
