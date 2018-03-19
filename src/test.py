import utilities

seqs = utilities.loadSequences("files/test/test.fasta", "files/test/test2.fasta")
for seq in seqs:
    print (seqs[seq].name)