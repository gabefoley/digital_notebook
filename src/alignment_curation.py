from Bio import AlignIO
import utilities
from collections import defaultdict


def get_seqs_at_index(alignment, col_index, return_gaps=True):
    """
    Given a column index get the sequences at that column that display a particular presence or absence of a gap
    :param alignment: The alignment to check
    :param col_index: The column indexes we're interested in
    :param gap: Whether we should return sequences that have a gap here
    :return: A list of sequences that display the particular presence or absence of a gap
    """
    seqs = []
    total_seqs = [x for x in range(len(alignment))]
    for align_index in range(len(alignment)):
        # If this column isn't a gap add it to insertion_columns list
        if alignment[align_index].seq[col_index] == "-":  #Add in all the sequences with gaps
            seqs.append(align_index)

    # If we want the sequences without gaps, just take the difference
    if not return_gaps:
        seqs = set(total_seqs) - set(seqs)
    return seqs


def check_sequential():
    pass


def check_position():
    pass

# for alignment in handle:
#     for idx in indexes:
#         for seq in alignment:
#             print (seq.seq[idx])


def get_deletion_candidates(alignment, min_length = 0, ignore_edges=True, rank=True, single=False) :
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
    for seq in alignment:
        print (seq.seq[4])
        if deletion in seq.seq:
            print ('candidate')
        # Order sequences with multiple deletions of minimum length first

        # Then order sequences by length of deletion

def get_insertion_candidates(alignment, min_length=0, ignore_edges=True, rank=True, single=False, accepted_percent = 0):
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


def get_indexes_of_deletions(alignment, accepted_percent):
    indexes = []
    # print(alignment)
    for idx in range(len(alignment[0])):
        col = []
        for seq in alignment:
            col += seq.seq[idx]
        col_deletions = (col.count("-"))
        if col_deletions > 0:
            percent_deletions = col_deletions / (len(col))
            # print(percent_deletions)
            if 1 - percent_deletions <= accepted_percent:
                indexes.append(idx)
    return indexes



def get_candidate_deletions(alignment, indexes, min_length):
    # print (indexes)
    candidates = {}
    positions = {}
    for i in range(len(indexes)):
        # print (i)
        if i+min_length <= len(indexes):
            window = indexes[i:i+min_length]
            idx = 0
            # print ('idx is now ', idx)
            while idx + 1 < len(window):

                # print ("window[idx] + 1 = ", window[idx] + 1)
                # print ("window[idx+1] = ", window[idx+1])
                # Check to see if the next column is a neighbouring column
                if window[idx] + 1 != window[idx + 1]:
                    idx = len(window)
                else:

                    # Get the list of alignments that have
                    insertion_cols = get_seqs_at_index(alignment, window[idx])
                    # print (insertion_cols)

                    if idx == 0:
                        prev_insertion_cols = intersection = insertion_cols
                        intersection = insertion_cols

                    else:
                        intersection = [val for val in prev_insertion_cols if val in insertion_cols]
                    if not intersection:
                        idx = len(window)

                # If we made it through the sliding window then it is a candidate insertion
                if idx +1 == len(window) - 1:
                        for seq in intersection:
                            candidate_index = 0
                            # print ("Looking at seq ", alignment[seq].name)
                            if alignment[seq].name not in candidates:
                                candidates[alignment[seq].name] = []
                                positions[alignment[seq].name] = []
                            else:
                                candidate_index = len(candidates[alignment[seq].name])
                            # candidates[alignment[seq].name][candidate_index].append(list(window))
                            first_index = window[0]
                            last_index = window[1]

                            if first_index not in positions[alignment[seq].name]:
                                # candidates[alignment[seq].name].append([])
                                if len(candidates[alignment[seq].name]) <= candidate_index:
                                    candidates[alignment[seq].name].append([])
                                candidates[alignment[seq].name][candidate_index].append(first_index)
                                positions[alignment[seq].name].append(first_index)
                                while first_index >= 0:
                                    if alignment[seq][first_index - 1] == "-":
                                        candidates[alignment[seq].name][candidate_index].append(first_index - 1)
                                        positions[alignment[seq].name].append(first_index - 1)
                                        first_index -= 1
                                    else:
                                        break
                            if last_index not in positions[alignment[seq].name]:
                                # candidates[alignment[seq].name].append([])
                                if len(candidates[alignment[seq].name]) <= candidate_index:
                                    candidates[alignment[seq].name].append([])
                                candidates[alignment[seq].name][candidate_index].append(last_index)
                                positions[alignment[seq].name].append(last_index)
                                # print ("last index")
                                # print (last_index)
                                # print ("len alignment + 1")
                                # print (len(alignment) +1)
                                while last_index <= len(alignment[0]) + 1 :
                                    # print ("alignment here is ", alignment[seq][last_index + 1])
                                    if alignment[seq][last_index + 1] == "-":
                                        candidates[alignment[seq].name][candidate_index].append(last_index + 1)
                                        positions[alignment[seq].name].append(last_index +1)
                                        # print ('candidates is ', candidates)
                                        last_index += 1
                                    else:
                                        break
                            # print (candidates)
                        #
                        #
                        # print ('candidates are', candidates)
                        # print("found window is", window)
                        # print("column is ", intersection)


                    # Expand the window in both possible directions

                        idx +=1
    print (candidates)
    return candidates

def get_seq_with_longest_deletion(candidate_deletions):
        lengths = defaultdict(list)
        for seq, positions in candidate_deletions.items():

            # Get the maximum length of a deletion for a single sequence
            max_len = len(max(positions, key=len))

            # If the maximum length is longer or equal to the current longest, add it to the dictionary
            if not lengths or max_len >= max(lengths, key=int):
                lengths[(max_len)].append(seq)

        #If there is only one candidate with the longest deletion, return it
        longest_candidates = lengths[(max(lengths, key=int))]
        print (longest_candidates)

        if len(longest_candidates) == 1:
            return longest_candidates[0]
        else:
            get_seq_with_most_deletions({k: candidate_deletions[k] for k in longest_candidates})


def get_seq_with_most_deletions(candidate_deletions):
    print (candidate_deletions)



alignment = utilities.load_alignment("../files/test/alignment_test.aln", "fasta")
accepted_percent = 0.9

indexes = get_indexes_of_deletions(alignment, accepted_percent)
print (alignment)
# print ("Indexes")
# print(indexes)

candidate_deletions = get_candidate_deletions(alignment, indexes, 2)

candidate_sequence = get_seq_with_longest_deletion(candidate_deletions)

print (candidate_sequence)

# print (candidate_deletions)
#
# print ("Indexes")
# print(indexes)