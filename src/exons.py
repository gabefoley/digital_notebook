from GenomicRecord import GenomicRecord
from Bio import SeqIO, Entrez
from collections import defaultdict
from dictsearch.search import iterate_dictionary
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
from reportlab.lib.units import cm
import fasta
import re
import math
import pickle
from urllib.request import urlopen
import utilities


Entrez.email = "gabriel.foley@uqconnect.edu.au"

# small = SeqIO.to_dict(SeqIO.parse("files/exons/small.fasta", "fasta"))
# homo_sapiens_broken = SeqIO.to_dict(SeqIO.parse("files/exons/homo_sapiens_broken.fasta", "fasta"))
# uniprot = SeqIO.to_dict(SeqIO.parse("files/exons/uniprot.fasta", "fasta"))
#
# latest_2U1 = SeqIO.to_dict(SeqIO.parse("/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and bacteria/2U1_55_percent.fasta", "fasta"))
# latest_2U1_cutdown = SeqIO.to_dict(SeqIO.parse("/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and bacteria/2U1_55_percent_cutdown.fasta", "fasta"))
# latest_2U1_2seqs = SeqIO.to_dict(SeqIO.parse("/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and bacteria/2U1_55_percent_2seqs.fasta", "fasta"))
# latest_2U1_3seqs = SeqIO.to_dict(SeqIO.parse("/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and bacteria/2U1_55_percent_3seqs.fasta", "fasta"))
# latest_2U1_3seqs = SeqIO.to_dict(SeqIO.parse("/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and bacteria/2U1_55_percent_3seqs.fasta", "fasta"))

def mapExons (records):
    """
    Take a set of protein sequences and map them back to their exon coordinates from their genomic location
    :param records: The set of protein sequences
    :return: The exon coordinates
    """

    genomic_records = {}

    for record in records:
        # print (record.id)
        try:
            protein_id = record.id

            if protein_id.startswith("tr"):
                protein_id = protein_id.split("|")[1]

                handle = urlopen("https://www.uniprot.org/uniprot/" + protein_id + ".xml")
                record = SeqIO.read(handle, "uniprot-xml")
                # print (record)
                for ref in record.dbxrefs:
                    if ref.startswith('EnsemblFungi'):
                        protein_id = (ref.split(":")[1])

                print (protein_id)

                # handle = Entrez.efetch(db="protein", id=gene_id, rettype="gb", retmode="xml")
                #
                # gene_record = Entrez.read(handle, validate=False)

                # print (gene_record)

            elif "|" in protein_id:
                protein_id = protein_id.split("|")[0]


            # Map the protein to a gene
            handle = Entrez.elink(dbfrom="protein", db='gene', id=protein_id)
            mapping = Entrez.read(handle, validate=False)

            print (mapping)


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
                print (protein_id)
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

                            print (feature_names)

                            for pos2 in range (0, len(feature_names)):
                                if feature_names[pos2]['GBQualifier_name'] == 'coded_by':
                                    i = feature_names[pos2]['GBQualifier_value']

                            # print ('koin?')
                            # print (i)

                            if 'join' in i:
                                exons = (i.split('join(')[1].split(','))

                                for count, exon in enumerate(exons):
                                    exons[count] = exon.split('.1:')[1]

                                strand = "minus" if "complement" in i else "plus"


                                # print ('The protein with ID %s has %d exons and thier locations are %s' % ( protein_id, len(exons), exons))
                                # print ('And the lengths of the exons are %s' % (update_exon_lengths(exons)))
                                # print ()

                                genomic_record = GenomicRecord(protein_id = protein_id, gene_id = gene_id, exons = exons, strand=strand, calc_introns =True)
                                genomic_records[protein_id] = genomic_record

            if gene_id:
                print ('and i got here')
                # These are the bits I removed but could add back
                # handle = Entrez.efetch(db="protein", id=protein_id, rettype="gb", retmode="xml")
                # full_record = Entrez.read(handle, validate=False)

                # for term in full_record:
                #     refseq = (term['GBSeq_source-db'].split("accession")[1].strip())

                # Use the gene ID to get the full gene record
                handle = Entrez.efetch(db="gene", id=gene_id, rettype="gb", retmode="xml")

                gene_record = Entrez.read(handle, validate=False)

                print (gene_record)



                # Get the genome ID, and the start and end location of the gene
                for term in gene_record:
                    print ('here now')
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
                                                print ('here')
                                                print (k[pos3])

                                                # print ('Equal one and in')

                                                # print ('%s had a Primary Assembly labelled genome' % protein_id)

                                                if 'Seq-loc_int' in k[pos3]['Gene-commentary_seqs'][0]:
                                                    genome_id = k[pos3]['Gene-commentary_accession']
                                                    seq_from = \
                                                        k[pos3]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval'][
                                                            'Seq-interval_from']
                                                    seq_to = \
                                                    k[pos3]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval'][
                                                        'Seq-interval_to']
                                                    strand_string = str(k[pos3]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_strand'])
                                                    result = re.search('attributes={\'value\': \'(.*)\'}\)}', strand_string)
                                                    strand = result.group(1)

                                            elif len(k) == 1 and 'Primary Assembly' not in k[pos3][
                                                'Gene-commentary_heading']:

                                                # print ('nope here')

                                                # print ('Equal one and not in')

                                                # print ('%s had a Primary Assembly labelled genome' % protein_id)

                                                if 'Seq-loc_int' in k[pos3]['Gene-commentary_seqs'][0]:
                                                    genome_id = k[pos3]['Gene-commentary_accession']
                                                    seq_from = \
                                                        k[pos3]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval'][
                                                            'Seq-interval_from']
                                                    seq_to = \
                                                    k[pos3]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval'][
                                                        'Seq-interval_to']
                                                    strand_string = str(k[pos3]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_strand'])
                                                    result = re.search('attributes={\'value\': \'(.*)\'}\)}', strand_string)
                                                    strand = result.group(1)
                                            elif len(k) > 1 and 'Primary Assembly' in k[pos3]['Gene-commentary_heading']:


                                                # print ('Greater than one and in')
                                                # if 'Seq-loc_int' in k[pos3]['Gene-commentary_seqs'][0] and 'Primary Assembly' in k[pos3]['Gene-commentary_heading']:
                                                if 'Seq-loc_int' in k[pos3]['Gene-commentary_seqs'][0]:
                                                    genome_id = k[pos3]['Gene-commentary_accession']
                                                    seq_from = \
                                                        k[pos3]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval'][
                                                            'Seq-interval_from']
                                                    seq_to = \
                                                    k[pos3]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval'][
                                                        'Seq-interval_to']
                                                    strand_string = str(k[pos3]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_strand'])
                                                    result = re.search('attributes={\'value\': \'(.*)\'}\)}', strand_string)
                                                    strand = result.group(1)



                    # print("The original protein record is %s which has a gene ID %s and a RefSeq ID %s " % (protein_id, gene_id, refseq))
                    print ('searching on')
                    print (genome_id)
                    # print ()
                    handle = Entrez.efetch(db="nuccore", id=genome_id, rettype="gb", retmode="xml", seq_start=seq_from,
                                           seq_stop=seq_to)


                    genome_record = Entrez.read(handle, validate=False)

                    for term in genome_record:
                        print (term)
                        feature_table = term['GBSeq_feature-table']
                        for pos in range (0, len(feature_table)):

                            if 'GBFeature_key' in feature_table[pos].keys() and feature_table[pos]['GBFeature_key'] == 'CDS':
                                feature_names = feature_table[pos]['GBFeature_quals']
                                for pos2 in range (0, len(feature_names)):

                                    print (feature_names[pos2]['GBQualifier_value'] )
                                    # Only pull the specific record we're searching for
                                    if feature_names[pos2]['GBQualifier_value'] == record.id:
                                        print ('here')
                                        print(feature_names[pos2]['GBQualifier_name'])
                                        print(feature_names[pos2]['GBQualifier_value'])

                                        protein_id = feature_names[pos2]['GBQualifier_value']

                                        i = feature_table[pos]['GBFeature_location']
                                        exons = (i.split('join(')[1].split(','))
                                        strand = "minus" if "complement" in i else "plus"

                                        # print ('The protein with ID %s has %d exons and thier locations are %s' % ( protein_id, len(exons), exons))
                                        # print ('And the lengths of the exons are %s' % (update_exon_lengths(exons)))
                                        # print ()

                                        # print (protein_id)
                                        # # print (gene_id)
                                        # print (strand)

                                        print ('exons are')
                                        print (exons)

                                        genomic_record = GenomicRecord(protein_id = protein_id, gene_id = gene_id, refseq = "",
                                                                       genome_id = genome_id, genome_positions =[seq_from, seq_to],
                                                                       exons = exons, strand=strand, calc_introns =True)
                                        # genomic_record.update_feature_lengths()

                                        print ('stored exons are')
                                        print (genomic_record.exons)
                                        print (genomic_record.exon_lengths)

                                        genomic_records[record.id] = genomic_record
                                        print ('added')
                                        print (genomic_records)

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


def appendTagToSeqsGivenExonCount(*numbers, exon_counts, full_record, outpath):
    exon_records = []
    outfile = []

    for count in numbers:
        if count in exon_counts:
            for id in exon_counts[count]:
                if id in full_record:
                    full_record[id].id += "*TAG " + str(count) + "*"

    for record in full_record.values():
        outfile.append(record)

    fasta.write_fasta(outfile, outpath)

def writeOutSeqsGivenExonCount(*numbers, exon_counts, full_record, outpath):

    exon_records = []
    for count in numbers:
        exon_records += exon_counts[count]

    print ('lets write it out')
    outfile = fasta.map_list_to_records(exon_records, full_record)

    fasta.write_fasta(outfile, outpath)

# records = mapExons(latest_2U1)
# counts = get_feature_counts(records)
# print (counts)
# appendTagToSeqsGivenExonCount(3,4,5,6,7,8,9,11,12,15, exon_counts=counts, full_record=latest_2U1, outpath="results/exons/2U1_55percent_tagged.fasta" )
#

# writeOutSeqsGivenExonCount(4,5,6, exon_counts=counts, full_record=latest_2U1, outpath="results/exons/2U1_55_percent_456_exons.fasta")

# for count, seqs in counts.items():
#     print ("Count: %d " % (count))
#     for seq in seqs:
#         print ("The exon lengths are %s " % (records[seq].exon_lengths))
#         print("The original protein record is %s which has a gene ID %s and a RefSeq ID %s " % (records[seq].protein_id, records[seq].gene_id, records[seq].refseq))
#         # print('The genome ID is %s and the gene region is from position %s to position %s ' % (records[seq].genome_id,
#         #                                                                                records[seq].genome_positions[0],
#         #                                                                                    records[seq].genome_positions[1]))
#         print ('The exons are %s and the exon lengths are %s' % (records[seq].exons, records[seq].exon_lengths))

# drawGenome()


def saveGenomicRecords(records, filepath):

    genomic_records = mapExons(records)
    utilities.savePythonObject(genomic_records, filepath)

def openGenomicRecords(filepath):
    genomic_records = utilities.openPythonObject(filepath)
    return genomic_records


def mapExonBoundariesToAlignment(records, genomic_records):

    # exons = {'XP_019684690.2' : [178, 167, 161, 635, 489], 'XP_009883824.1' : [178, 167, 161, 705, 82], 'XP_014065111.1' : [181, 167, 161, 334, 318, 441] }
    # exons = {'XP_019684690.2' : [3, 3, 2, 1], 'XP_009883824.1' : [2, 2, 1], 'XP_014065111.1' : [4] }
    #


    # genomic_records = mapExons(records)



    cols = ["\033[1;35;4m", "\033[1;34;4m", "\033[1;32;4m", "\033[1;33;4m",  "\033[1;37;4m" ,  "\033[1;31;4m" ,  "\033[1;39;4m"  ]
    longest_header = 0
    for record in records:
        if len(record.id) > longest_header:
            longest_header = len(record.id)


        buildseq = ""
        newseq = ""

        counter = 0

        # print (genomic_records[record].strand)
        # print (genomic_records[record].exon_lengths)
        # print (genomic_records[record].exons)

        # skip_records = ['XP_010997363.1', 'XP_006162319.2', 'XP_005406339.2', 'XP_014050304.1'] # This was the skip records for 55% identity
        skip_records = ['XP_010997363.1', 'XP_021107396.1', 'XP_002607780.1', 'PIK57532.1', 'PIK38219.1', 'XP_002605102.1', 'XP_006825012.1', 'XP_014801266.1', 'XP_021107397.1', 'XP_006162319.2', 'PIO41157.1', 'XP_005406339.2', 'KYO44822.1', 'KTF79201.1', 'XP_014379039.1', 'XP_013865517.1',  'XP_006520454.1', 'ELW64418.1', 'XP_016094679.1', 'EPY81189.1', 'XP_016404253.1', 'XP_006768464.1', 'XP_013766542.1', 'XP_006123186.1', 'XP_014050304.1', 'ETE67196.1', 'KTF79021.1', 'OBS74422.1', 'XP_004671492.1', 'EMP40787.1', 'KFO25848.1','KKF09825.1'] # This was the skip records for 50% identity

        # skip_records = ['XP_010224405.1', 'XP_014050304.1', 'XP_006162319.2', 'XP_010997363.1', 'XP_005406339.2']
        # skip_records = ['XP_010224405.1', 'XP_006162319.2', 'XP_010997363.1', 'XP_005406339.2']

        # print (record.id)
        # print (genomic_records)

        if record.id in genomic_records and record.id not in skip_records :

            # print (genomic_records[record].exon_lengths)
            # print (genomic_records[record].strand)
            for exon in genomic_records[record.id].exon_lengths:
                # print (exon)
                # print (exon / 3)
                counter += (exon/3)
            # print (counter)

            for count, exon in enumerate(genomic_records[record.id].exon_lengths):
                # print (count, exon)
                # print (count,int(math.ceil(exon/3)))
                buildseq += str(count+1) * int(math.ceil((exon / 3)))

            # print (record)
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
                    ind +=1

            # print (newseq)

            # print (records[record].seq)
            seq_with_exons =  str(record.seq)
            # print (newseq)
            # print (len(genomic_records[record].exon_lengths))
            for exon in range (0, len(genomic_records[record.id].exon_lengths)):
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


# records_4seqs = SeqIO.to_dict(SeqIO.parse(
#     "/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and bacteria/2U1_55_percent_4seqs.aln",
#     "fasta"))
#
# latest_2U1 = SeqIO.to_dict(SeqIO.parse(
#     "/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and bacteria/2U1_latest.aln",
#     "fasta"))
#
# new_2U1 = list(SeqIO.parse(
#     "/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and bacteria/2U1_latest.aln",
#     "fasta"))
#
# # latest_2U1 = SeqIO.to_dict(SeqIO.parse(
# #     "/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and bacteria/problem.fasta",
# #     "fasta"))
#
# latest_2U1_3seqs = SeqIO.to_dict(SeqIO.parse(
#     "/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and bacteria/test.aln",
#     "fasta"))
#
# old_2U1 = list(SeqIO.parse("/Users/gabefoley/Dropbox/Code/Python Workspace/digital_notebook/files/2U 2R aligned.fasta", "fasta"))
#
# original_2U1 = list(SeqIO.parse("/Users/gabefoley/Dropbox/Code/Python Workspace/digital_notebook/files/raine_exons.fasta", "fasta"))
#
#
# bug = list(SeqIO.parse("/Users/gabefoley/Dropbox/Code/Python Workspace/digital_notebook/files/bug.fasta", "fasta"))
#
#
# leo = list(SeqIO.parse("/Users/gabefoley/Dropbox/PhD/Smaller Projects/Leo Errors/uniprot_single.fasta", "fasta"))



# genome_mapping = list(SeqIO.parse("/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and bacteria/180312_fifty_percent_identity/Genome_mapping/2U1_55percent_tagged.aln",
#     "fasta"))





# filepath = 'problem.obj'
genome_mapping = list(SeqIO.parse("/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and bacteria/180312_fifty_percent_identity/Genome_mapping/2U1_50percent.aln",
    "fasta"))
filepath = 'files/objects/genome_mapping.obj'
# saveGenomicRecords(genome_mapping, filepath)

genomic_records = openGenomicRecords(filepath)

# for x in genomic_records:
#     print (x)

mapExonBoundariesToAlignment(genome_mapping, genomic_records)

fifty_percent = SeqIO.to_dict(SeqIO.parse("/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and bacteria/180312_fifty_percent_identity/Genome_mapping/2U1_50_percent.fasta",
    "fasta"))

# counts = get_feature_counts(genomic_records)
# print (counts.keys())
# appendTagToSeqsGivenExonCount(3,4,5,6,7,8,9,10,11,12,13,15, exon_counts=counts, full_record=fifty_percent, outpath="results/exons/2U1_50percent_tagged.fasta" )
# writeOutSeqsGivenExonCount(3,4,5,6,7,8,9,10,11,12,13,15, exon_counts=counts, full_record=fifty_percent, outpath="results/exons/2U1_50_percent_tagged.fasta")

# errors = list(SeqIO.parse("/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and bacteria/180312_fifty_percent_identity/Genome_mapping/errors.aln",
#     "fasta"))
#
#
# filepath = 'files/objects/errors.obj'
#
# saveGenomicRecords(errors, filepath)
#
# genomic_records = openGenomicRecords(filepath)
#
# for x in genomic_records:
#     print (x)
#
# mapExonBoundariesToAlignment(errors, genomic_records)



