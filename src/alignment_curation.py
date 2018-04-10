from Bio import AlignIO
import utilities
from collections import defaultdict


def get_seqs_at_index(alignment, col_index, filter=None, return_gaps=True):
    """
    Given a column index get the sequences at that column that display a particular presence or absence of a gap
    :param alignment: The alignment to check
    :param col_index: The column indexes we're interested in
    :param filter: The list of sequences to only include
    :param return_gaps: Whether we should return sequences that have a gap here
    :return: A list of sequences that display the particular presence or absence of a gap
    """
    seqs = []
    # print ('FILTER IS ', filter)
    total_seqs = [x for x in range(len(alignment))]
    if filter:
        check_seqs = filter
    else:
        check_seqs = [x for x in range(len(alignment))]
    # print ('check seqs is ', check_seqs)
    for align_index in check_seqs:

        # If this column isn't a gap add it to insertion_columns list
        if alignment[align_index].seq[col_index] == "-":  # Add in all the sequences with gaps
            seqs.append(align_index)
        else:
            if align_index in seqs:
                seqs.pop(align_index)

    # If we want the sequences without gaps, just take the difference
    if not return_gaps:
        seqs = set(total_seqs) - set(seqs)
    return seqs


def check_sequential():
    pass


def check_position(alignment, seq, index):
    pass


# for alignment in handle:
#     for idx in indexes:
#         for seq in alignment:
#             print (seq.seq[idx])


# def get_deletion_candidates(alignment, min_length=0, ignore_edges=True, rank=True, single=False):
#     """
#     Return a list of the candidates to be removed based on presence of deletion
#     :param alignment: The alignment file to search
#     :param min_length: The minimum length of the deletion
#     :param ignore_edges: Whether to ignore deletions occuring at the N or C terminus
#     :param rank: Whether to rank sequences by most / longest deletions
#     :param single: Whether to just return a single candidate
#     :return: List of candidates / candidate to be removed base on presence of a deletion
#     """
#     deletion = "-" * min_length
#     print(deletion)
#     for seq in alignment:
#         if deletion in seq.seq:
#             print('candidate')
#             # Order sequences with multiple deletions of minimum length first
#
#             # Then order sequences by length of deletion


def get_insertion_candidates(alignment, min_length=0, ignore_edges=True, rank=True, single=False, accepted_percent=0):
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


def get_indexes_of_deletions(alignment, accepted_percent, filter=None):
    indexes = []
    for idx in range(len(alignment[0])):
        col = alignment[:,idx]

        if filter:

            passed_filter = False
            for filter_seq in filter:
                if col[filter_seq] == "-":
                    passed_filter = True
                    break

        if not filter or passed_filter:

            col_deletions = (col.count("-"))
            if col_deletions > 0:
                percent_deletions = col_deletions  / (len(col) )

                # Check if it meets the minimum percent deletion for a column and if the gap is present in one of the sequences we're filtering on
                # print (col)
                # print ("percent_deletions ",  1 - percent_deletions)
                # print ("accepted_percent ", accepted_percent)
                if 1 - percent_deletions >= accepted_percent:
                    indexes.append(idx)
    return indexes


def get_candidate_deletions(alignment, indexes, min_length, filter=None, internal_only=True):
    # print (indexes)
    candidates = {}
    positions = {}
    for i in range(len(indexes)):
        if i + min_length <= len(indexes):
            window = indexes[i:i + min_length]
            idx = 0
            # print('idx is now ', idx)
            # print(len(window))

            # If the user just wants a window size of one
            if len(window) == 1:
                insertion_cols = get_seqs_at_index(alignment, window[0], filter)
                intersection = [val for val in insertion_cols]


            else:
                while idx + 1  <= len(window):

                    # print ("window[idx] + 1 = ", window[idx] + 1)
                    # print ("window[idx+1] = ", window[idx+1])
                    # Check to see if the next column is a neighbouring column
                    if idx + 1 < len(window) and window[idx] + 1 != window[idx + 1]:
                        idx = len(window)
                    else:

                        # Get the list of alignments that have
                        insertion_cols = get_seqs_at_index(alignment, window[idx], filter)
                        # print ('idx is ', window[idx])
                        # print ("insertion cols", insertion_cols)

                        if idx == 0:
                            prev_insertion_cols = intersection = insertion_cols
                            intersection = insertion_cols

                        else:
                            intersection = [val for val in prev_insertion_cols if val in insertion_cols]
                        if not intersection:
                            idx = len(window)

                    # print ('intersection', intersection)
                    idx += 1


            # idx = 0
            # If we made it through the sliding window then it is a candidate insertion
            if idx  == len(window):
                for seq in intersection:
                    internal = True

                    candidate_index = 0
                    if (alignment[seq].name == "cat"):
                        print ("cat")
                    ("Looking at seq ", alignment[seq].name)
                    if alignment[seq].name not in candidates:
                        candidates[alignment[seq].name] = []
                        positions[alignment[seq].name] = []
                    else:
                        candidate_index = len(candidates[alignment[seq].name])
                    # candidates[alignment[seq].name][candidate_index].append(list(window))
                    first_index = window[0]
                    last_index = window[-1]

                    block_filled = False

                    if first_index not in positions[alignment[seq].name]:
                        # candidates[alignment[seq].name].append([])
                        if len(candidates[alignment[seq].name]) <= candidate_index:
                            candidates[alignment[seq].name].append([])
                        first_region = window[0:-1]

                        candidates[alignment[seq].name][candidate_index].extend(first_region)
                        positions[alignment[seq].name].extend(first_region)
                        while first_index >= 0:
                            # If we only want internal deletions, check that we haven't extended to the N-terminal
                            if internal_only:
                                if first_index == 0:
                                    candidates[alignment[seq].name].pop(candidate_index)
                                    internal = False
                                    break



                            if alignment[seq][first_index - 1] == "-":
                                candidates[alignment[seq].name][candidate_index].insert(0, first_index - 1)
                                positions[alignment[seq].name].append(first_index - 1)
                                first_index -= 1
                            else:
                                break
                    if last_index not in positions[alignment[seq].name] and (not internal_only or internal):
                                # candidates[alignment[seq].name].append([])
                                if len(candidates[alignment[seq].name]) <= candidate_index:
                                    candidates[alignment[seq].name].append([])
                                candidates[alignment[seq].name][candidate_index].append(last_index)
                                positions[alignment[seq].name].append(last_index)
                                # print ("last index")
                                # print (last_index)
                                # print ("len alignment ")
                                # print (len(alignment[0]))
                                while last_index + 1  < len(alignment[0]):

                                    # print ("alignment here is ", alignment[seq][last_index + 1])
                                    if alignment[seq][last_index + 1] == "-":
                                        candidates[alignment[seq].name][candidate_index].append(last_index + 1)
                                        positions[alignment[seq].name].append(last_index + 1)
                                        # print ('candidates is ', candidates)
                                        last_index += 1
                                    else:
                                        break
                                # if last_index + 1 == len(alignment[0]):
                                #     if alignment[seq][last_index] == "-":
                                #         candidates[alignment[seq].name][candidate_index].append(last_index + 1)
                                #         positions[alignment[seq].name].append(last_index + 1)




    remove = {}

    # Check for minimum length and minimum coverage
    for seq, locations in candidates.items():
        for index, location in enumerate(locations):
            if (len(location) < min_length):
                remove[seq] = index

    for seq, pos in remove.items():
        candidates[seq].pop(pos)

    # Check if this made any of the candidate lists empty
    remove = [k for k in candidates if len(candidates[k]) == 0]
    for k in remove: del candidates[k]

    # print ('REMOVER', remove)

   # remove_seq = [] for pos in candidates.values():
   #      if (len(pos)) == 0:


    # print(candidates)

    return candidates


def get_seq_with_longest_deletion(candidate_deletions, final_check = True):

    lengths = defaultdict(list)
    if candidate_deletions:
        for seq, positions in candidate_deletions.items():


            # Get the maximum length of a deletion for a single sequence
            max_len = len(max(positions, key=len))

            # If the maximum length is longer or equal to the current longest, add it to the dictionary
            if not lengths or max_len >= max(lengths, key=int):
                lengths[(max_len)].append(seq)

        print(lengths)

        if lengths:
            longest_candidates = lengths[(max(lengths, key=int))]

            # Check if we want to use longest deletion as the final check or if we only have one candidate
            if final_check or len(longest_candidates) == 1:
                return longest_candidates
            else:
                candidate = get_seq_with_most_deletions({k: candidate_deletions[k] for k in longest_candidates})
                return candidate
        else:
            return None
    else:
        return None


def get_seq_with_most_deletions(candidate_deletions, final_check=True):
    print ('checking for most deletions')
    nums = defaultdict(list)
    for seq, positions in candidate_deletions.items():
        # print(seq)
        # print(len(positions))
        max_num = len(positions)
        if not nums or max_num >= max(nums, key=int):
            nums[(max_num)].append(seq)

    if nums:
        # print (nums)
        highest_num_candidates = nums[(max(nums, key=int))]
        # print (highest_num_candidates)
        # print (highest_num_candidates[0])

        if final_check or len(highest_num_candidates) == 1:
            return highest_num_candidates
        else:
            get_seq_with_longest_deletion({k: candidate_deletions[k] for k in highest_num_candidates})
    else:
        return None


def filter_alignment_by_deletion_length(alignment,length, internal_only=True):
    """
    Filter an alignment so that we return a list of sequences that contain at least one deletion of the minimum length
    :param alignment: The alignment to check
    :param length: The minimum length
    :return: A list of sequences that contain a deletion of minimum length
    """
    sequences = []
    for seq_index in range(len(alignment)):
        if internal_only:
            seq = alignment[seq_index].seq.strip("-")

        else:
            seq = alignment[seq_index]

        if "-" * length in seq:
            sequences.append(seq_index)
    return sequences

