from Bio import AlignIO, SeqIO

records = SeqIO.parse("files/example_expected.fasta", "fasta")

def checkLongIndels(record, indelLength):
    checkLongSubstring
    if longestSubString > indelLength:
        return True
    else:
        return False


for record in records:
    checkLongIndels(record.seq, 5)
    print (alignment)

