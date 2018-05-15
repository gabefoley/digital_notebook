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
import csv

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

        print('Search id is ')
        print(search_id)

        # If it isn't an NCBI sequence, lets try and map it to the UniProt database
        if record.annotations["Database"] != "NCBI":

            try:
                handle = urlopen("https://www.uniprot.org/uniprot/" + search_id + ".xml")
                soup = Soup(handle, "lxml")
                search_id = soup.find('property', type="protein sequence ID")["value"]

            except HTTPError as error:
                pass
                print ("Couldn't access a Uniprot record for " + search_id)

        if protein_record:
            # print ('Here is the NCBI id')
            # print (search_id)

            # Check if this is an ID that maps to the Ensembel Fungi Genome database
            if (search_id[:4] == "FOXG"):
                handle = urlopen("http://www.ensemblgenomes.org/id/" + search_id)
                soup = Soup(handle, "lxml")

                transcript_info = soup.find('div', attrs={'class': 'lhs'}, text="About this transcript")

                for sibling in transcript_info.next_siblings:
                    if sibling != "\n":
                        exon_num = int(sibling.text.split("This transcript has ")[1].split(" ")[0])

                exons = []

                for i in range(0, exon_num):
                    exons.append("id.1:1..1")

                # print (exons)

                for count, exon in enumerate(exons):
                    exons[count] = re.split("[.]\d:", exon)[1]

                strand = "plus"

                genomic_record = GenomicRecord(protein_id=search_id, gene_id="", exons=exons,
                                               strand=strand, calc_introns=False)

                # print ("Exon count is", genomic_record.exon_count)
                genomic_records[record.id] = genomic_record

            else:
                handle = Entrez.efetch(db="protein", id=search_id, rettype="gb", retmode="xml")
                soup = Soup(handle, "lxml")
                coded_by = soup.find("gbqualifier_name", text="coded_by")

                for sibling in coded_by.next_siblings:
                    if sibling != "\n":
                        exon_text = sibling.text
                        # print (exon_text)

                # Check that the record isn't coded by mRNA
                if exon_text[0:2] != "XM":

                    if 'join' in exon_text:
                        exons = (exon_text.split('join(')[1].split(','))

                        for count, exon in enumerate(exons):
                            exons[count] = re.split("[.]\d:", exon)[1]



                    else:
                        exons = [re.split("[.]\d:", exon_text)[1]]

                    strand = "minus" if "complement" in exon_text else "plus"

                    genomic_record = GenomicRecord(protein_id=search_id, gene_id="", exons=exons,
                                                   strand=strand, calc_introns=True)

                    print("Exon count is", genomic_record.exon_count)

                    genomic_records[record.id] = genomic_record

                # The record was coded by mRNA, so we need to try and map the Protein record to a Gene
                else:

            #     # Map the protein to a gene
            #
                    handle = Entrez.elink(dbfrom="protein", db='gene', id=search_id, rettype='xml')
                    print ('This was an mRNA record')

                    try:
                        mapping = Entrez.read(handle, validate=False)
                    except:
                        mapping = None
                        gene_id = None

                    # Retrieve the gene ID
                    if mapping:
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
                                            exons[count] = re.split("[.]\d:", exon)[1]



                                    else:
                                        print('no join')
                                        print(i)
                                        exons = [re.split("[.]\d:", i)[1]]


                                    print ('find me an exon', exons)

                                    strand = "minus" if "complement" in i else "plus"

                                    genomic_record = GenomicRecord(protein_id=search_id, gene_id=gene_id, exons=exons,
                                                                   strand=strand, calc_introns=True)
                                    genomic_records[record.id] = genomic_record

                    if gene_id:
                        # print('We could map the NCBI protein ID to a gene automatically')
                        # print(gene_id)

                        # Use the gene ID to get the full gene record
                        handle = Entrez.efetch(db="gene", id=gene_id, rettype="gb", retmode="xml")
                        # print (handle.read())
                        # gene_record = Entrez.read(handle, validate=False)
                        # print (gene_record)
                        soup = Soup(handle, "lxml")
                        # print (soup.prettify())
                        seq_from = soup.find("seq-interval_from")

                        seq_to = soup.find("seq-interval_to")
                        genome_id_list = soup.find_all('gene-commentary_accession')

                        for id in genome_id_list:
                            id_text = id.getText()
                            # print (id_text)
                            if id_text[0:2] == 'NW':
                                genome_id = id_text


                        # print("The original protein record is %s which has a gene ID %s and a RefSeq ID %s "
                        # % (search_id, gene_id, refseq))
                        # print('searching on')
                        # print(genome_id)
                        # print ()

                        if (genome_id):
                            handle = Entrez.efetch(db="nuccore", id=genome_id, rettype="gb", retmode="xml",
                                                   seq_start=seq_from,
                                                   seq_stop=seq_to)
                            soup = Soup(handle, "lxml")
                            # print(soup.prettify())
                            # print (gene_id)
                            gb_features = soup.find_all('gbfeature_key', text="CDS")

                            for feature in gb_features:
                                # print (feature)
                                for qualifier in feature.find_all_next('gbqualifier_value'):
                                    gene_check = qualifier.getText().split(":")
                                    if gene_check[0] == "GeneID" and gene_check[1] == gene_id:
                                        # print ('found it')
                                        for location in feature.find_all_next('gbfeature_location'):
                                            exon_location = (location.getText())

                            if 'join' in exon_location:
                                # print('joint is here')
                                # print (exon_location)
                                exons = (exon_location.split('join(')[1].split(','))

                                # print (exons)

                                # for count, exon in enumerate(exons):
                                #     exons[count] = re.split("[.]\d:", exon)[1]



                            else:
                                # print('no join')
                                # print(exon_location)
                                exons = [re.split("[.]\d:", exon_location)[1]]

                            # print('find me an exon', exons)

                            strand = "minus" if "complement" in exon_location else "plus"

                            genomic_record = GenomicRecord(protein_id=search_id, gene_id=gene_id,
                                                           exons=exons,
                                                           strand=strand, calc_introns=True)

                            print ("Exon count is ", genomic_record.exon_count)
                            genomic_records[record.id] = genomic_record




                        else:
                            print ("Couldn't find a genome ID")


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

    # print(genomic_records)
    utilities.save_python_object(genomic_records, filepath)


def open_genomic_records(filepath):
    genomic_records = utilities.open_python_object(filepath)
    return genomic_records


def get_exon_counts(records, genomic_records, filter_records=None, exclude=True):
    exon_counts = {}
    if not filter_records:
        filter_records = []

    if type(records) == dict:
        records = [record for record in records.values()]

    for record in records:

        if record.id in genomic_records:
            if (exclude and record.id not in filter_records) or (not exclude and record.id in filter_records):
                exon_counts[record.id] = len(genomic_records[record.id].exons)
    return exon_counts


def write_exon_counts_to_csv(records, filepath, totals=True):
    with open(filepath, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=' ',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        csvwriter.writerow("Hello| Goodby")
        # csvwriter.writerow(['Sequence ID']['Exon Count'])
        # for records, exons in records.items():
        #     csvwriter.writerow([records][exons])
        # if totals:
        #     csvwriter.writerow["-"]["Total"]


def write_exon_totals_to_csv(records, filepath):
    with open(filepath, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=' ',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        csvwriter.writerow("Hello| Goodby")


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

