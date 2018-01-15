from GenomicRecord import GenomicRecord
from Bio import SeqIO, Entrez
from collections import defaultdict
from dictsearch.search import iterate_dictionary

Entrez.email = "gabriel.foley@uqconnect.edu.au"

small = SeqIO.to_dict(SeqIO.parse("files/exons/small.fasta", "fasta"))
homo_sapiens_broken = SeqIO.to_dict(SeqIO.parse("files/exons/homo_sapiens_broken.fasta", "fasta"))
uniprot = SeqIO.to_dict(SeqIO.parse("files/exons/uniprot.fasta", "fasta"))

def mapExons (records):
    """
    Take a set of protein sequences and map them back to their exon coordinates from their genomic location
    :param records: The set of protein sequences
    :return: The exon coordinates
    """

    genomic_records = {}

    for record in records:
        try:
            protein_id = records[record].id

            if "|" in protein_id:
                protein_id = protein_id.split("|")[0]

            print (protein_id)

            # Map the protein to a gene
            handle = Entrez.elink(dbfrom="protein", db='gene', id=protein_id)
            mapping = Entrez.read(handle, validate=False)


            # Retrieve the gene ID
            for term in mapping:
                if (term['LinkSetDb']):
                    if (term['LinkSetDb'][0]['Link'][0]['Id']):
                        gene_id = (term['LinkSetDb'][0]['Link'][0]['Id'])
                else:
                    gene_id = None
                    break
            if gene_id == None:
                print ('was none')
                handle = Entrez.efetch(db="protein", id=protein_id, rettype="gb", retmode="xml")


                genome_record = Entrez.read(handle, validate=False)
                # print (genome_record)

                for term in genome_record:
                    feature_table = term['GBSeq_feature-table']
                    # print (feature_table[0].keys())
                    for pos in range (0, len(feature_table)):
                        if 'GBFeature_key' in feature_table[pos].keys() and feature_table[pos]['GBFeature_key'] == 'CDS':
                            # print (feature_table[pos])
                            # i = feature_table[pos]['GBFeature_location']
                            feature_names = feature_table[pos]['GBFeature_quals']

                            # print (feature_names)

                            for pos2 in range (0, len(feature_names)):
                                if feature_names[pos2]['GBQualifier_name'] == 'coded_by':
                                    i = feature_names[pos2]['GBQualifier_value']

                            # print ('koin?')
                            # print (i)

                            exons = (i.split('join(')[1].split(','))


                            for count, exon in enumerate(exons):
                                exons[count] = exon.split('.1:')[1]


                            # print ('The protein with ID %s has %d exons and thier locations are %s' % ( protein_id, len(exons), exons))
                            # print ('And the lengths of the exons are %s' % (update_exon_lengths(exons)))
                            # print ()

                            genomic_record = GenomicRecord(protein_id = protein_id, gene_id = gene_id, exons = exons, calc_introns =True)
                            genomic_records[protein_id] = genomic_record

            else:
                handle = Entrez.efetch(db="protein", id=protein_id, rettype="gb", retmode="xml")
                full_record = Entrez.read(handle, validate=False)

                for term in full_record:
                    refseq = (term['GBSeq_source-db'].split("accession")[1].strip())

                # Use the gene ID to get the full gene record
                handle = Entrez.efetch(db="gene", id=gene_id, rettype="gb", retmode="xml")

                gene_record = Entrez.read(handle, validate=False)



                # Get the genome ID, and the start and end location of the gene
                for term in gene_record:
                    i = term['Entrezgene_comments']
                    for pos in range(0, len(i)):
                        if 'Gene-commentary_comment' in i[pos].keys():
                            j = i[pos]['Gene-commentary_comment']
                            for pos2 in range(0, len(j)):
                                if 'Gene-commentary_products' in j[pos2].keys():
                                    k = j[pos2]['Gene-commentary_products']

                                    for pos3 in range(0, len(k)):
                                        if 'Gene-commentary_accession' in k[pos3].keys() and 'Gene-commentary_seqs' in k[
                                            pos3].keys() and 'Gene-commentary_heading' in k[pos3].keys():

                                            if len(k) == 1 and 'Primary Assembly' in k[pos3]['Gene-commentary_heading']:
                                                print ('%s had a Primary Assembly labelled genome' % protein_id)

                                                if 'Seq-loc_int' in k[pos3]['Gene-commentary_seqs'][0]:
                                                    genome_id = k[pos3]['Gene-commentary_accession']
                                                    seq_from = \
                                                        k[pos3]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval'][
                                                            'Seq-interval_from']
                                                    seq_to = \
                                                    k[pos3]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval'][
                                                        'Seq-interval_to']

                                            elif len(k) == 1 and 'Primary Assembly' not in k[pos3][
                                                'Gene-commentary_heading']:

                                                print ('%s had a Primary Assembly labelled genome' % protein_id)

                                                if 'Seq-loc_int' in k[pos3]['Gene-commentary_seqs'][0]:
                                                    genome_id = k[pos3]['Gene-commentary_accession']
                                                    seq_from = \
                                                        k[pos3]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval'][
                                                            'Seq-interval_from']
                                                    seq_to = \
                                                        k[pos3]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval'][
                                                            'Seq-interval_to']
                                            elif len(k) > 1 and 'Primary Assembly' in k[pos3]['Gene-commentary_heading']:
                                                # if 'Seq-loc_int' in k[pos3]['Gene-commentary_seqs'][0] and 'Primary Assembly' in k[pos3]['Gene-commentary_heading']:
                                                if 'Seq-loc_int' in k[pos3]['Gene-commentary_seqs'][0]:
                                                    print('%s had a Primary Assembly labelled genome' % protein_id)

                                                    genome_id = k[pos3]['Gene-commentary_accession']
                                                    seq_from = \
                                                        k[pos3]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval'][
                                                            'Seq-interval_from']
                                                    seq_to = \
                                                        k[pos3]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval'][
                                                            'Seq-interval_to']




                    # print("The original protein record is %s which has a gene ID %s and a RefSeq ID %s " % (protein_id, gene_id, refseq))

                    handle = Entrez.efetch(db="nuccore", id=genome_id, rettype="gb", retmode="xml", seq_start=seq_from,
                                           seq_stop=seq_to)


                    genome_record = Entrez.read(handle, validate=False)

                    for term in genome_record:
                        feature_table = term['GBSeq_feature-table']
                        for pos in range (0, len(feature_table)):
                            if 'GBFeature_key' in feature_table[pos].keys() and feature_table[pos]['GBFeature_key'] == 'CDS':
                                i = feature_table[pos]['GBFeature_location']
                                feature_names = feature_table[pos]['GBFeature_quals']

                                for pos2 in range (0, len(feature_names)):
                                    if feature_names[pos2]['GBQualifier_name'] == 'protein_id':
                                        protein_id = feature_names[pos2]['GBQualifier_value']

                                exons = (i.split('join(')[1].split(','))
                                # print ('The protein with ID %s has %d exons and thier locations are %s' % ( protein_id, len(exons), exons))
                                # print ('And the lengths of the exons are %s' % (update_exon_lengths(exons)))
                                # print ()

                                genomic_record = GenomicRecord(protein_id = protein_id, gene_id = gene_id, refseq = refseq,
                                                               genome_id = genome_id, genome_positions =[seq_from, seq_to],
                                                               exons = exons, calc_introns =True)
                                genomic_records[protein_id] = genomic_record

                                # print('The genome ID is %s' % ( genome_id))
                                # print('The exons are %s and the exon lengths are %s' % (exons, genomic_record.exon_lengths))
                                #
                                #


                                # exon_records.append(exon_dict)
        except Exception:
            continue
    return genomic_records

def get_feature_counts(records):
    """
    Takes an exon records dictionary and returns a dictionary with exon count as the key and a list of the speices that
    have that number exons as the value

    :param exon_records: The dictionary containing all of the exon records
    :return: Dictionary with the exon counts as the key
    """

    exon_counts = defaultdict(list)

    for protein_id, genomic_record in records.items():
        exon_counts[genomic_record.exon_count].append(protein_id)
    return exon_counts


def map_exon_count_to_tree(exons, tree):
    print ('done')

# records = mapExons(homo_sapiens_broken)
# records = mapExons(small)
# records = mapExons(uniprot)
#
#
# counts = get_feature_counts(records)
#
# for count, seqs in counts.items():
#     print ("Count: %d " % (count))
#     for seq in seqs:
#         print ("The exon lengths are %s " % (records[seq].exon_lengths))
#         print("The original protein record is %s which has a gene ID %s and a RefSeq ID %s " % (records[seq].protein_id, records[seq].gene_id, records[seq].refseq))
#         # print('The genome ID is %s and the gene region is from position %s to position %s ' % (records[seq].genome_id,
#         #                                                                                records[seq].genome_positions[0],
#         #                                                                                    records[seq].genome_positions[1]))
#         print ('The exons are %s and the exon lengths are %s' % (records[seq].exons, records[seq].exon_lengths))