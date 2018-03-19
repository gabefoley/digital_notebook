import fasta
import utilities


seqs = utilities.loadSequences("files/exons/uniprot.fasta")

fasta.print_record_overview(seqs)