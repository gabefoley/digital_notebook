import utilities
import fasta

seqs = utilities.loadSequences("files/exons/uniprot.fasta")

seqs = fasta.add_species_from_header(seqs)

for seq in seqs.values():
    print(seq.name)
    print(seq.annotations["Species"])
