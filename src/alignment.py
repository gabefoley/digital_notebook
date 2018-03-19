from Bio.Align.Applications import MafftCommandline, ClustalOmegaCommandline
from Bio import AlignIO
import io


def align_with_mafft(filepath):
    """
    Align a file with the given filepath using MAFFT
    :param filepath: The file to align
    :return: The MAFFT alignment
    """
    mafft_cline = MafftCommandline(input=filepath)
    stdout, stderr = mafft_cline()
    align = AlignIO.read(io.StringIO(stdout), "fasta")

    return align


def align_with_clustal_omega(filepath):
    """
    Align a file with the given filepath using Clustal Omega
    :param filepath: The file to align
    :return: The Clustal Omega alignment
    """
    clustal_omega_cline = ClustalOmegaCommandline(infile=filepath)

    stdout, stderr = clustal_omega_cline()
    align = AlignIO.read(io.StringIO(stdout), "fasta")
    return align


def write_alignment(file, filepath, filetype):
    """
    Write an alignment to file
    
    :param file: The alignment
    :param filepath: The filepath to write to
    :param filetype: The type of file to create
    """
    AlignIO.write(file, filepath, filetype)


def get_percent_identity(seq1, seq2, count_gaps=False):
    """
    Calculate the percent identity between two aligned sequences
    :param seq1: First sequence
    :param seq2: Second sequence
    :param count_gaps: Whether we should count a gap as a mismatch or no
    :return:
    """
    matches = sum(aa1 == aa2 for aa1, aa2 in zip(seq1, seq2))

    # Set the length based on whether we want identity to count gaps or not
    length = len(seq1) if count_gaps else min(len(seq1.replace("-", "")), len(seq2.replace("-", "")))

    pct_identity = 100.0 * matches / length

    return pct_identity
