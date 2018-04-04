import utilities
import exons

filepath = "../files/uniprot/exons_15.obj"
# uniprot_seqs = utilities.load_sequences("../files/uniprot/uniprot_small.fasta", split_char="|")

uniprot_seqs = utilities.load_sequences("../files/uniprot/uniprot_15.fasta", split_char="|")

print ("There are this many uniprot seqs in the original file - %d" % (len(uniprot_seqs)))

exons.save_genomic_records(uniprot_seqs, filepath)
genomic_records = exons.open_genomic_records(filepath)

exons.map_exon_boundaries_to_alignment(uniprot_seqs, genomic_records)