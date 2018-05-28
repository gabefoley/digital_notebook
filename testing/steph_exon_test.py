import utilities
import fasta
import exons
import glob
from Bio import SeqIO
import sys

working_dir = "/Users/gabefoley/Dropbox/PhD/Smaller Projects/Stephlina_exons"

files = glob.glob(working_dir + "/Al6onAl5.fasta")
sequences = utilities.load_sequences(working_dir + "/Alignment_1.fasta")
sequences = utilities.load_sequences(working_dir + "/test.fasta")


files = glob.glob(working_dir + "/no_automatic_exons.fasta")


save_files = True
sys.setrecursionlimit(50000)
total_counts = {}

for file in files:
    name = file.rsplit("/")[-1].split(".")[0]
    filepath = working_dir + "/objects/" + name + ".obj"
    skip_path = working_dir + "/output/" + name + "_skipped.csv"


    print("Working on ", name)
    if save_files:
        print("Saving to ", filepath)
    alignment_file = utilities.load_sequences(file, split_char="|")
    if save_files:
        exons.save_genomic_records(alignment_file, filepath, skipped_records_path=skip_path)

    exon_dict = exons.open_genomic_records(filepath)

    exons.write_exon_counts_to_csv(exon_dict, working_dir + "/output/" + name + "_exons.csv")
    #
    # exon_array = exons.get_exon_array(exon_dict)
    # total_counts[name] = exon_array


exons.write_exon_totals_to_csv(total_counts, working_dir + "/output/total_exons.csv")



exons.map_exon_boundaries_to_alignment(sequences, exon_dict, filter_records=["X5KCU7", "X5KNV9", "B7ZP10", "H5ZVH5",
                                                                             "B0BMP5", "A9JTU3", "Q7Z449", "G9L0D0", "K9IST0", "Q0IIF9", "A8KBF5"])