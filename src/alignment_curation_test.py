import alignment_curation
import utilities

alignment = utilities.load_alignment("../files/test/small_test.fasta", "fasta")


# alignment = utilities.load_alignment("../files/test/large.aln", "fasta")
alignment = utilities.load_alignment("../files/test/2U1_50_percent_motif.aln", "fasta")
accepted_percent = 0.2
min_length = 30
internal_only = True

# print (alignment)
# print (" ")


filters = alignment_curation.filter_alignment_by_deletion_length(alignment, min_length, internal_only)
#
# print ("Filters are ", filters)
#
indexes = alignment_curation.get_indexes_of_deletions(alignment, accepted_percent, filters)

# print ("Indexes are ", indexes)
#
candidate_deletions = alignment_curation.get_candidate_deletions(alignment, indexes, min_length, filters, internal_only)
# print ("Candidate deletions are ", candidate_deletions)
#
candidate_sequence = alignment_curation.get_seq_with_longest_deletion(candidate_deletions, final_check=False)

print ("The candidate sequence is")
print(candidate_sequence)
['XP_023051618.1', 'XP_023372804.1*TAG', 'XP_016807465.1*TAG', 'XP_017708555.1*TAG', 'XP_017202860.1*TAG', 'XP_018880964.1*TAG']
