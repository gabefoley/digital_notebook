from Bio import SeqIO, Entrez
from collections import defaultdict
Entrez.email = "gabriel.foley@uqconnect.edu.au"
homo_sapiens = SeqIO.to_dict(SeqIO.parse("files/exons/homo_sapiens.fasta", "fasta"))
strongylocentrotus_purpuratus = SeqIO.to_dict(SeqIO.parse("files/exons/strongylocentrotus_purpuratus.fasta", "fasta"))
all_seqs = SeqIO.to_dict(SeqIO.parse("files/exons/all.fasta", "fasta"))
agaricus_bisporus = SeqIO.to_dict(SeqIO.parse("files/exons/Agaricus_bisporus_stripped_removed.fasta", "fasta"))
filtered_2U1 = SeqIO.to_dict(SeqIO.parse("files/exons/2U1_filtered_records.fasta", "fasta"))
filtered_pruned_2U1 = SeqIO.to_dict(SeqIO.parse("files/exons/2U1_BLAST_filtered_records_X_removed_C_2.fasta", "fasta"))


class GenomicRecord(object):

    protein_id = ""
    nucleotides = ""
    amino_acids = ""
    exons = []
    exon_lengths = []
    intron_lengths = []

    def __init__(self, protein_id, gene_id = "", refseq="", genome_id = "",  genome_positions = [], exons = [], introns = [], calc_introns=False):
            self.protein_id = protein_id
            self.gene_id = gene_id
            self.refseq = refseq
            self.genome_id = genome_id
            self.genome_positions = genome_positions
            self.nucleotides = ""
            self.amino_acids = ""
            self.exons = exons
            self.introns = introns
            self.exon_count = len(exons)
            self.exon_lengths = []
            self.intron_lengths = []

            self.update_feature_lengths(calc_introns)




    def get_introns_and_exons(self):
        regions = []
        for region in zip(self.introns, self.exons):
            regions.append(region)
        return regions

    def update_feature_lengths(self, calc_introns=False):
        # print (self.exons)
        for count, exon in enumerate(self.exons):

            exon_start = int(exon.split("..")[0].replace(")", "").replace("<", "").replace(">", ""))
            exon_end = int(exon.split("..")[1].replace(")", "").replace("<", "").replace(">", ""))

            if calc_introns:
                if count != 0 and count != len(self.exons):
                    curr_exon = exon_start
                    self.intron_lengths.append(curr_exon - prev_exon)

                prev_exon = exon_end

            self.exon_lengths.append(exon_end - exon_start)

