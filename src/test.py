import utilities

seqs = utilities.load_sequences("files/test/test.fasta", "files/test/test2.fasta")
for seq in seqs:
    print (seqs[seq].name)