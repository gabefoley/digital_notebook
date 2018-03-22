import fasta
import utilities


seqs = utilities.load_sequences("files/exons/uniprot.fasta")

fasta.print_record_overview(seqs)