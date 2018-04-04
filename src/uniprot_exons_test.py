import utilities
import exons

filepath = "../files/uniprot/exons_small.obj"
uniprot_seqs = utilities.load_sequences("../files/uniprot/uniprot_small.fasta", split_char="|")

# uniprot_seqs = utilities.load_sequences("../files/uniprot/uniprot_large.fasta")

for seq in uniprot_seqs.values():
    print (seq.name)
    print (seq.id)
    print (" ")

print ("There are this many uniprot seqs in the original file - %d" % (len(uniprot_seqs)))

exons.saveGenomicRecords(uniprot_seqs, filepath)
genomic_records = exons.openGenomicRecords(filepath)


for record, value in genomic_records.items():
    print (record)
    print (value)

exons.mapExonBoundariesToAlignment(uniprot_seqs, genomic_records)