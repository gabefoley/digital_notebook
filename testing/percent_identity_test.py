import alignment
import utilities

alignment_file = utilities.load_alignment("/Users/gabefoley/Dropbox/Code/Python Workspace/digital_notebook/testing/alignment.fasta")

percent_identities_array = alignment.get_percent_identity_of_alignment(alignment_file)

print("Percent identity is ", percent_identities_array.mean())
print("Standard deviation is ", percent_identities_array.std())
