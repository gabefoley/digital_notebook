import utilities
import fasta

seqs = utilities.loadSequences("files/exons/uniprot.fasta")

seqs = fasta.addSpeciesFromHeader(seqs)

for seq in seqs.values():
    print (seq.name)
    print (seq.annotations["Species"])

