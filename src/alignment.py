from Bio import AlignIO, SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
import io
import re

filepath = "files/2U1_all_candidates_PSI_BLAST.fasta"


def alignWithMAFFT(filepath):
    mafft_cline = MafftCommandline(input=filepath)
    stdout, stderr = mafft_cline()
    align = AlignIO.read(io.StringIO(stdout), "fasta")

    return align

def writeAlignment(file, filepath, format):
    AlignIO.write(file, filepath, format)

def checkLongIndels(record, indelLength):
    longestIndel = len(max(re.compile("(-+-)*").findall(str(record))))
    if longestIndel > indelLength:
        return True
    else:
        return False



# mafftAlign = alignWithMAFFT(filepath)


# print (mafftAlign)
# for line in mafftAlign:
#     print (line)
#     print (checkLongIndels(line.seq, 5))

# def getPercentIdentity(seq1, seq2):
#