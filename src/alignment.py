from Bio import AlignIO, SeqIO
from Bio.Align.Applications import MafftCommandline, ClustalOmegaCommandline
from Bio import AlignIO
import io
import re

filepath = "files/priapulus_caudatus.fasta"


def alignWithMAFFT(filepath):
    mafft_cline = MafftCommandline(input=filepath)
    stdout, stderr = mafft_cline()
    align = AlignIO.read(io.StringIO(stdout), "fasta")

    return align


def align_with_clustal_omega(filepath):
    clustal_omega_cline = ClustalOmegaCommandline(infile=filepath)

    stdout, stderr = clustal_omega_cline()
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

def getPercentIdentity(seq1, seq2, count_gaps=False):
    matches = sum(aa1 == aa2 for aa1, aa2 in zip(seq1, seq2))
    # Set the length based on whether we want identity to count gaps or not

    length = len(seq1) if count_gaps else min(len(seq1.replace("-", "")), len(seq2.replace("-", "")))

    pct_identity = 100.0 * matches / length
    return pct_identity

# clustal_align = align_with_clustal_omega(filepath)

# for line in clustal_align:
#     print (line)