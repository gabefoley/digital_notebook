from GenomicRecord import GenomicRecord
from Bio import SeqIO, Entrez
from collections import defaultdict
import fasta
import re
import math
from urllib.request import urlopen
import utilities
from bs4 import BeautifulSoup as Soup
from urllib.error import HTTPError

Entrez.email = "gabriel.foley@uqconnect.edu.au"


def map_exons(records):
    """
    Take a set of protein sequences and map them back to their exon coordinates from their genomic location
    :param records: The set of protein sequences
    :return: The exon coordinates
    """
    genomic_records = {}

    for record in records:
        # try:
        search_id = record.id
        protein_record = True

        # If it isn't an NCBI sequence, lets try and map it to the UniProt database
        if record.annotations["Database"] != "NCBI":
            print ('search id is ')
            print (search_id)
            try:
                handle = urlopen("https://www.uniprot.org/uniprot/" + search_id + ".xml")
                soup = Soup(handle, "lxml")
                search_id = soup.find('property', type="protein sequence ID")["value"]
            except HTTPError as error:
                print ("Couldn't access a Uniprot record for " + search_id)

        if protein_record:
            # Map the protein to a gene
            handle = Entrez.elink(dbfrom="protein", db='gene', id=search_id)
            mapping = Entrez.read(handle, validate=False)

            print('This is the protein id')
            print(search_id)

            # print ('This is the mapping')
            # print(mapping)

            # Retrieve the gene ID
            for term in mapping:
                if term['LinkSetDb']:
                    if term['LinkSetDb'][0]['Link'][0]['Id']:
                        gene_id = term['LinkSetDb'][0]['Link'][0]['Id']
                else:
                    gene_id = None
                    break

            if gene_id is None:
                print("Couldn't map the NCBI protein ID to a gene automatically")
                print(search_id)
                handle = Entrez.efetch(db="protein", id=search_id, rettype="gb", retmode="xml")

                genome_record = Entrez.read(handle, validate=False)
                # print (genome_record)

                for term in genome_record:
                    feature_table = term['GBSeq_feature-table']
                    # print (feature_table[0].keys())
                    for pos in range(0, len(feature_table)):
                        if 'GBFeature_key' in feature_table[pos].keys() and \
                                        feature_table[pos]['GBFeature_key'] == 'CDS':
                            # print (feature_table[pos])
                            # i = feature_table[pos]['GBFeature_location']
                            feature_names = feature_table[pos]['GBFeature_quals']

                            print(feature_names)

                            for pos2 in range(0, len(feature_names)):
                                if feature_names[pos2]['GBQualifier_name'] == 'coded_by':
                                    i = feature_names[pos2]['GBQualifier_value']

                            print(i)

                            if 'join' in i:
                                print('joint is here')
                                print(i)
                                exons = (i.split('join(')[1].split(','))

                                for count, exon in enumerate(exons):
                                    exons[count] = exon.split('.1:')[1]

                            else:
                                print('no join')
                                print(i)
                                exons = [i.split('.1:')[1]]

                            print(exons)

                            strand = "minus" if "complement" in i else "plus"

                            genomic_record = GenomicRecord(protein_id=search_id, gene_id=gene_id, exons=exons,
                                                           strand=strand, calc_introns=True)
                            genomic_records[record.id] = genomic_record

            if gene_id:
                print('We could map the NCBI protein ID to a gene automatically')
                print(gene_id)

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
                                        if 'Gene-commentary_accession' in k[
                                            pos3].keys() and 'Gene-commentary_seqs' in \
                                                k[
                                                    pos3].keys() and 'Gene-commentary_heading' in k[pos3].keys():

                                            if len(k) == 1 and 'Primary Assembly' in \
                                                    k[pos3]['Gene-commentary_heading']:
                                                if 'Seq-loc_int' in k[pos3]['Gene-commentary_seqs'][0]:
                                                    genome_id = k[pos3]['Gene-commentary_accession']
                                                    seq_from = \
                                                        k[pos3]['Gene-commentary_seqs'][0]['Seq-loc_int'][
                                                            'Seq-interval'][
                                                            'Seq-interval_from']
                                                    seq_to = \
                                                        k[pos3]['Gene-commentary_seqs'][0]['Seq-loc_int'][
                                                            'Seq-interval'][
                                                            'Seq-interval_to']
                                                    strand_string = str(
                                                        k[pos3]['Gene-commentary_seqs'][0]['Seq-loc_int'][
                                                            'Seq-interval']['Seq-interval_strand'])
                                                    result = re.search('attributes={\'value\': \'(.*)\'}\)}',
                                                                       strand_string)
                                                    strand = result.group(1)

                                            elif len(k) == 1 and 'Primary Assembly' not in \
                                                    k[pos3]['Gene-commentary_heading']:

                                                if 'Seq-loc_int' in k[pos3]['Gene-commentary_seqs'][0]:
                                                    genome_id = k[pos3]['Gene-commentary_accession']
                                                    seq_from = \
                                                        k[pos3]['Gene-commentary_seqs'][0]['Seq-loc_int'][
                                                            'Seq-interval'][
                                                            'Seq-interval_from']
                                                    seq_to = \
                                                        k[pos3]['Gene-commentary_seqs'][0]['Seq-loc_int'][
                                                            'Seq-interval'][
                                                            'Seq-interval_to']
                                                    strand_string = str(
                                                        k[pos3]['Gene-commentary_seqs'][0]['Seq-loc_int'][
                                                            'Seq-interval']['Seq-interval_strand'])
                                                    result = re.search('attributes={\'value\': \'(.*)\'}\)}',
                                                                       strand_string)
                                                    strand = result.group(1)
                                            elif len(k) > 1 and 'Primary Assembly' in \
                                                    k[pos3]['Gene-commentary_heading']:

                                                # print ('Greater than one and in')
                                                # if 'Seq-loc_int' in k[pos3]['Gene-commentary_seqs'][0] and
                                                # 'Primary Assembly' in k[pos3]['Gene-commentary_heading']:
                                                if 'Seq-loc_int' in k[pos3]['Gene-commentary_seqs'][0]:
                                                    genome_id = k[pos3]['Gene-commentary_accession']
                                                    seq_from = \
                                                        k[pos3]['Gene-commentary_seqs'][0]['Seq-loc_int'][
                                                            'Seq-interval'][
                                                            'Seq-interval_from']
                                                    seq_to = \
                                                        k[pos3]['Gene-commentary_seqs'][0]['Seq-loc_int'][
                                                            'Seq-interval'][
                                                            'Seq-interval_to']
                                                    strand_string = str(
                                                        k[pos3]['Gene-commentary_seqs'][0]['Seq-loc_int'][
                                                            'Seq-interval']['Seq-interval_strand'])
                                                    result = re.search('attributes={\'value\': \'(.*)\'}\)}',
                                                                       strand_string)
                                                    strand = result.group(1)

                    # print("The original protein record is %s which has a gene ID %s and a RefSeq ID %s "
                    # % (search_id, gene_id, refseq))
                    print('searching on')
                    print(genome_id)
                    # print ()
                    handle = Entrez.efetch(db="nuccore", id=genome_id, rettype="gb", retmode="xml",
                                           seq_start=seq_from,
                                           seq_stop=seq_to)

                    genome_record = Entrez.read(handle, validate=False)

                    for term in genome_record:
                        # print(term)
                        feature_table = term['GBSeq_feature-table']
                        for pos in range(0, len(feature_table)):

                            if 'GBFeature_key' in feature_table[pos].keys() and \
                                            feature_table[pos]['GBFeature_key'] == 'CDS':
                                feature_names = feature_table[pos]['GBFeature_quals']
                                # print (feature_names)
                                for pos2 in range(0, len(feature_names)):

                                    print("here is the name", feature_names[pos2]['GBQualifier_value'])
                                    if ":" in feature_names[pos2]['GBQualifier_value']:
                                        print (feature_names[pos2]['GBQualifier_value'].split(":")[1])
                                        print (gene_id)
                                    # Only pull the specific record we're searching for
                                        if feature_names[pos2]['GBQualifier_value'].split(":")[1] == gene_id:
                                            # print('here')
                                            # print(feature_names[pos2]['GBQualifier_name'])
                                            # print(feature_names[pos2]['GBQualifier_value'])
                                            print ("Here is the record id")
                                            print (feature_names[pos2])

                                            search_id = feature_names[pos2]['GBQualifier_value']

                                            i = feature_table[pos]['GBFeature_location']

                                            print(i)
                                            if "," in i:
                                                exons = (i.split('join(')[1].split(','))
                                            else:
                                                exons = i
                                            strand = "minus" if "complement" in i else "plus"

                                            # print ('The protein with ID %s has %d exons and thier locations are %s' %
                                            # ( search_id, len(exons), exons))
                                            # print ('And the lengths of the exons are %s' % (update_exon_lengths(exons)))
                                            # print ()

                                            # print (search_id)
                                            # # print (gene_id)
                                            # print (strand)

                                            # print('exons are')
                                            # print(exons)

                                            genomic_record = GenomicRecord(protein_id=search_id, gene_id=gene_id,
                                                                           refseq="",
                                                                           genome_id=genome_id,
                                                                           genome_positions=[seq_from, seq_to],
                                                                           exons=exons, strand=strand,
                                                                           calc_introns=True)
                                            # genomic_record.update_feature_lengths()
                                            #
                                            print('stored exons are')
                                            print(genomic_record.exons)
                                            # print(genomic_record.exon_lengths)

                                            genomic_records[record.id] = genomic_record
                                            # print('added')
                                            # print(genomic_records)

                                            # print('The genome ID is %s' % ( genome_id))
                                            # print('The exons are %s and the exon lengths are %s' %
                                            # (exons, genomic_record.exon_lengths))
                                            #
                                            #


                                            # exon_records.append(exon_dict)

                                            # else:
                                            #     print ("This is the gene id")
                                            #     print (search_id)
                                            #     # search_id = "KL657831.1"
                                            #
                                            #     handle = Entrez.efetch(db="nucleotide", id=search_id, rettype="gb",
                                            # retmode="xml")
                                            #     mapping = Entrez.read(handle, validate=False)
                                            #     # print('This is the mapping')
                                            #     # print(mapping)
                                            #     # print(mapping[0].keys())
                                            #     #
                                            #     # for k, v in mapping[0].items():
                                            #     #     print (k)
                                            #     #     print (" ")
                                            #     #     print (v)
                                            #
                                            #     coding_seq_count = 0
                                            #
                                            #
                                            #
                                            #     if "GBSeq_feature-table" in mapping[0].keys():
                                            #         feature_table = mapping[0]["GBSeq_feature-table"]
                                            #
                                            #         for pos in range(0, len(feature_table)):
                                            #             if feature_table[pos]["GBFeature_key"] == 'CDS':
                                            #                 # print ('yep')
                                            #                 feature_names = feature_table[pos]
                                            #                 print ("HERES DA LOCATION")
                                            #                 print (feature_names["GBFeature_location"])
                                            #                 coding_seq_count +=1
                                            #
                                            #     if coding_seq_count > 1:
                                            #         print ("**** WARNING **** Multiple coding sequences found")

        # except Exception:
        #     continue
    return genomic_records


def get_feature_counts(records):
    """
    Takes an exon records dictionary and returns a dictionary with exon count as the key and a list of the speices that
    have that number exons as the value

    :param records: The dictionary containing all of the exon records
    :return: Dictionary with the exon counts as the key
    """

    exon_counts = defaultdict(list)

    for protein_id, genomic_record in records.items():
        exon_counts[genomic_record.exon_count].append(protein_id)
    return exon_counts


def map_exon_count_to_tree(exons, tree):
    pass


def append_tag_to_seqs_given_exon_count(*numbers, exon_counts, full_record, outpath):
    outfile = []

    for count in numbers:
        if count in exon_counts:
            for seq_id in exon_counts[count]:
                if seq_id in full_record:
                    full_record[seq_id].id += "*TAG " + str(count) + "*"

    for record in full_record.values():
        outfile.append(record)

    fasta.write_fasta(outfile, outpath)


def write_out_seqs_given_exon_count(*numbers, exon_counts, full_record, outpath):
    exon_records = []
    for count in numbers:
        exon_records += exon_counts[count]

    outfile = fasta.map_list_to_records(exon_records, full_record)

    fasta.write_fasta(outfile, outpath)


def save_genomic_records(records, filepath):
    # If records is a dictionary convert it into a list of the records
    if type(records) == dict:
        records = [record for record in records.values()]

    genomic_records = map_exons(records)

    print(genomic_records)
    utilities.save_python_object(genomic_records, filepath)


def open_genomic_records(filepath):
    genomic_records = utilities.open_python_object(filepath)
    return genomic_records


def print_exon_counts(records, genomic_records, filter_records=None, exclude=True):
    if not filter_records:
        filter_records = []

    if type(records) == dict:
        records = [record for record in records.values()]

    for record in records:

        if record.id in genomic_records:
            print (record.id)
            if (exclude and record.id not in filter_records) or (not exclude and record.id in filter_records):
                #

                # print (genomic_records[record].exon_lengths)
                # print (genomic_records[record].strand)
                print( len(genomic_records[record.id].exons))


def map_exon_boundaries_to_alignment(records, genomic_records, filter_records=None, exclude=True):

    if not filter_records:
        filter_records = []

    # If records is a dictionary convert it into a list of the records
    if type(records) == dict:
        records = [record for record in records.values()]

    # exons = {'XP_019684690.2' : [178, 167, 161, 635, 489], 'XP_009883824.1' : [178, 167, 161, 705, 82],
    # 'XP_014065111.1' : [181, 167, 161, 334, 318, 441] }
    # exons = {'XP_019684690.2' : [3, 3, 2, 1], 'XP_009883824.1' : [2, 2, 1], 'XP_014065111.1' : [4] }
    #

    # genomic_records = map_exons(records)

    cols = ["\033[1;35;4m", "\033[1;34;4m", "\033[1;32;4m", "\033[1;33;4m", "\033[1;37;4m", "\033[1;31;4m",
            "\033[1;39;4m"]
    longest_header = 0
    for record in records:
        # print (record.id)
        if len(record.id) > longest_header:
            longest_header = len(record.id)

        buildseq = ""
        newseq = ""

        counter = 0

        # print (genomic_records[record].strand)
        # print (genomic_records[record].exon_lengths)
        # print (genomic_records[record].exons)




        if record.id in genomic_records:
            # print (record.id)
            if (exclude and record.id not in filter_records) or (not exclude and record.id in filter_records):
                print ('here')
                #

                # print (genomic_records[record].exon_lengths)
                # print (genomic_records[record].strand)
                for exon in genomic_records[record.id].exon_lengths:
                    # print (exon)
                    # print (exon / 3)
                    counter += (exon / 3)
                # print (counter)

                for count, exon in enumerate(genomic_records[record.id].exon_lengths):
                    # print (count, exon)
                    # print (count,int(math.ceil(exon/3)))
                    buildseq += str(count + 1) * int(math.ceil((exon / 3)))

                # print (buildseq)
                # print (record.seq)
                print(record.id)

                ind = 0
                for pos in record.seq:
                    if pos == "-":
                        newseq += "-"
                    elif pos != "-":
                        # print (ind)
                        newseq += buildseq[ind]
                        ind += 1

                # print (newseq)

                # print (records[record].seq)
                seq_with_exons = str(record.seq)
                # print (newseq)
                # print (len(genomic_records[record].exon_lengths))
                for exon in range(0, len(genomic_records[record.id].exon_lengths)):
                    # print ('ind', ind)
                    # print ('exon', exon + 1)
                    # print ('done')
                    ind = (newseq.index(str(exon + 1)))

                    # print (len(seq_with_exons))
                    seq_with_exons = seq_with_exons[:ind + exon] + "*" + seq_with_exons[ind + exon:]
                    # seq_with_exons =  str(exon) + "*" + seq_with_exons[:ind + 2 * exon]

                for exon in range(0, len(genomic_records[record.id].exon_lengths)):
                    seq_with_exons = seq_with_exons.replace("*", cols[exon], 1)
                    # print ('hereeeee')
                    # print (seq_with_exons)

                print(cols[-1] + '{message: <{fill}}'.format(message=record.id, fill=longest_header), seq_with_exons)

                # else:
                # print (cols[-1] + '{message: <{fill}}'.format(message=record.id, fill=longest_header), record.seq)
