import utilities
import fasta
import exons
import glob
from Bio import SeqIO

working_dir = "/Users/gabefoley/Dropbox/Code/Python Workspace/digital_notebook/notebooks/Latest exon testing/All_Full_aligned_cleaned_together"

files = glob.glob(working_dir + "/*.fasta")

save_files = True

for file in files:
    name = file.rsplit("/")[-1].split(".")[0]
    filepath = working_dir + "/objects/" + name + ".obj"

    print("Working on ", name)
    print("Saving to ", filepath)
    alignment_file = utilities.load_sequences(file, split_char="|")
    if save_files:
        exons.save_genomic_records(alignment_file, filepath)

    exon_dict = exons.open_genomic_records(filepath)
    exons.write_exon_counts_to_csv(exon_dict, working_dir + "/output/" + name + "_exons.csv")
    exons.write_exon_totals_to_csv(exon_dict, working_dir + "/output/" + name + "_total_exons.csv")

#     percent_identities_array = alignment.get_percent_identity_of_alignment(alignment_file)

#     print("Percent identity is ", percent_identities_array.mean())
#     print("Standard deviation is ", percent_identities_array.std())