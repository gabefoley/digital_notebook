import utilities
import fasta
import exons
import glob
from Bio import SeqIO
import sys
working_dir = "/Users/gabefoley/Dropbox/Code/Python Workspace/digital_notebook/notebooks/Latest exon testing/All_Full_aligned_cleaned_together"
working_dir = "/Users/gabefoley/Dropbox/All_Full_aligned_cleaned_together"

files = glob.glob(working_dir + "/*.fasta")

save_files = True
sys.setrecursionlimit(50000)
total_counts = {}

for file in files:
    name = file.rsplit("/")[-1].split(".")[0]
    filepath = working_dir + "/objects/" + name + ".obj"
    skip_path = working_dir + "/output/" + name + "_skipped.csv"


    print("Working on ", name)
    alignment_file = utilities.load_sequences(file, split_char="|")
    if save_files:
        print("Saving to ", filepath)

        exons.save_genomic_records(alignment_file, filepath, skipped_records_path=skip_path)

    exon_dict = exons.open_genomic_records(filepath)

    print ('exon dict is')
    print (exon_dict)

    if exon_dict:
        exons.write_exon_counts_to_csv(exon_dict, working_dir + "/output/" + name + "_exons.csv")

        exon_array = exons.get_exon_array(exon_dict)
        total_counts[name] = exon_array


exons.write_exon_totals_to_csv(total_counts, working_dir + "/output/total_exons.csv")

#     percent_identities_array = alignment.get_percent_identity_of_alignment(alignment_file)

#     print("Percent identity is ", percent_identities_array.mean())
#     print("Standard deviation is ", percent_identities_array.std())