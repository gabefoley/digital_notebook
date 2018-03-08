from Bio import AlignIO, SeqIO, Entrez
import re


handle = list(AlignIO.parse(
    "/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and bacteria/2U1_55_percent_4seqs.aln",
    "fasta"))

accepted_percent = 0
min_length = 16
indexes = []

def get_insertion_columns_at_index(alignment_list, col_index):
    for alignment in alignment_list:
        insertion_columns = []
        for align_index in range(len(alignment)):
            # If this column isn't a gap add it to insertion_columns list
            if alignment[align_index].seq[col_index] != "-":
                insertion_columns.append(align_index)
    return insertion_columns

for alignment in handle:
    for idx in range(len(alignment[0])):
        col = []
        for seq in alignment:
            col += seq.seq[idx]
        col_deletions = (col.count("-"))
        if col_deletions > 0:
            percent_deletions = col_deletions / (len(col) - 1)
            if 1 - percent_deletions <= accepted_percent:
                indexes.append(idx)
print (indexes)

n = 16
for i in range(len(indexes)):
    if i+n <= len(indexes):
        window = indexes[i:i+n]
        idx = 0
        while (idx + 1 < len(window)):
            # Check to see if the next column is a neighbouring column
            if (window[idx] + 1 != window[idx + 1]):
                idx = len(window)
            else:

                # Get the list of alignments that have
                insertion_cols = get_alignment_columns_at_index(handle, window[idx])

                if (idx == 0):
                    prev_insertion_cols = intersection = insertion_cols
                    intersection = insertion_cols

                else:
                    intersection = [val for val in prev_insertion_cols if val in insertion_cols]
                if not intersection:
                    idx = len(window)



            # If we made it through the sliding window then it is a candidate insertion
            if idx +1 == len(window) - 1:
                    print ("found window is", window)
                    print ("column is ", intersection)

                # Expand the window in both possible directions

            idx +=1



def check_sequential():
    pass

def check_position():
    pass

# for alignment in handle:
#     for idx in indexes:
#         for seq in alignment:
#             print (seq.seq[idx])

def get_deletion_candidates(handle, min_length = 0, ignore_edges=True, rank=True, single=False) :
    """
    Return a list of the candidates to be removed based on presence of deletion
    :param alignment: The alignment file to search
    :param min_length: The minimum length of the deletion
    :param ignore_edges: Whether to ignore deletions occuring at the N or C terminus
    :param rank: Whether to rank sequences by most / longest deletions
    :param single: Whether to just return a single candidate
    :return: List of candidates / candidate to be removed base on presence of a deletion
    """
    deletion = "-" * min_length
    print (deletion)
    for alignment in handle:
        for seq in alignment:
            print (seq.seq[4])
            if deletion in seq.seq:
                print ('candidate')
        # Order sequences with multiple deletions of minimum length first

        # Then order sequences by length of deletion

def get_insertion_candidates(handle, min_length=0, ignore_edges=True, rank=True, single=False, accepted_percent = 0):
    """
    Return a list of the candidates to be removed based on presence of an insertion
    :param alignment: The alignment file to search
    :param min_length: The minimum length of the insertion
    :param ignore_edges: Whether to ignore insertions occuring at the N or C terminus
    :param rank: Whether to rank sequences by most / longest insertions
    :param single: Whether to just pull a single candidate
    :param accepted_percent: The accepted percent of other sequences having character state at the column of the insertion
    :return: List of candidates / candidate to be removed base on presence of an insertion
    """


get_deletion_candidates(handle, 20)

