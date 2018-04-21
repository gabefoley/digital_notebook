import utilities
import fasta
import exons
from Bio import SeqIO

working_dir ="/Users/gabefoley/Dropbox/Code/Python Workspace/digital_notebook/notebooks/Latest exon testing/"

# Load in the file and define a filepath that we'll save the object onto

aao = utilities.load_sequences(working_dir + "AAO_Full_211-seq.fasta", split_char="|")

gox = utilities.load_sequences(working_dir + "GOX_Full_88-seq.fasta", split_char="|")
# gox = utilities.load_sequences(working_dir + "gox_error.fasta", split_char="|")


# aao = utilities.load_sequences(working_dir + "aao_small.fasta", split_char="|")

# gox = utilities.load_sequences(working_dir +"gox_small.fasta", split_char="|")

fail = utilities.load_sequences(working_dir + "fail_one.fasta", split_char="|")
fail = utilities.load_sequences(working_dir + "ncbi_looks_like_uniprot.fasta", split_char="|")


aao_filepath = working_dir + 'files/objects/AAO.obj'
gox_filepath = working_dir + 'files/objects/GOX.obj'
fail_filepath = working_dir + "files/objects/fail.obj"

# for seq in aao.values():
#     print (seq.id)

# Save the files
# exons.save_genomic_records(aao, aao_filepath)
exons.save_genomic_records(gox, gox_filepath)
# exons.save_genomic_records(fail, fail_filepath)

# Open the file
# aao_records = exons.open_genomic_records(aao_filepath)
# print (len(aao_records))
# # exons.map_exon_boundaries_to_alignment(cyp2U1_exons, genomic_records, filter_records=["XP_010997363.1", "XP_006162319.2", "XP_004671492.1", "XP_005406339.2"], exclude=True)
# exons.print_exon_counts(aao, aao_records, exclude=True)

# Open the file
gox_records = exons.open_genomic_records(gox_filepath)
print (gox_records)
# exons.map_exon_boundaries_to_alignment(cyp2U1_exons, genomic_records, filter_records=["XP_010997363.1", "XP_006162319.2", "XP_004671492.1", "XP_005406339.2"], exclude=True)
exons.print_exon_counts(gox, gox_records, exclude=True)

# Open the file
# fail_records = exons.open_genomic_records(fail_filepath)
# print (fail_records)
# # exons.map_exon_boundaries_to_alignment(cyp2U1_exons, genomic_records, filter_records=["XP_010997363.1", "XP_006162319.2", "XP_004671492.1", "XP_005406339.2"], exclude=True)
# exons.print_exon_counts(fail, fail_records, exclude=True)