import utilities
import fasta
import re

string = "> find this..."

search = re.search('>(.*)\.\.\.', string)
print (search.group(1).strip())


seqs = utilities.loadSequences("/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and bacteria/180312_fifty_percent_identity/Alignment_exclusion/Comparison_to_old_alignment/2U1_50_percent_motif_9.fasta")
cd_hit = open('/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and bacteria/180312_fifty_percent_identity/Alignment_exclusion/Comparison_to_old_alignment/2U1_50_percent_motif_9_cdhit100.fasta.clstr')

fasta.checkCDHitOutput(seqs, cd_hit)