from Bio import SeqIO, Entrez
from collections import defaultdict
# import plotly
# plotly.tools.set_credentials_file(username='gabefoley', api_key='xS8qT0kIbIKDWt0BalOd')
# import plotly.plotly as py
# import plotly.graph_objs as go
# from plotly.graph_objs import *
import re

out_records =[]

# Builds up a dictionary with species as keys and lists of record ids that are from that species as values
# def build_species_count(*exclude_on, records, length=0 ):
#     species_count = defaultdict(list)
#     for record in records.values():
#         if (len(str(record.seq)) > length and not any(exclude in record.description for exclude in exclude_on)):
#             if "[" in record.description:
#                 name = record.description.split(">")[0].split("[")[1].split("]")[0]
#             else:
#                 name = "_".join(record.description.split("|")[1].split("_")[0:2])
#             species_count[name].append(record.id)
#     return species_count



def subset_records(*header_terms, records, length=0, mode="exclude"):
    """

    :param header_terms:
    :param mode:
    :return:
    """
    new_records = {}
    for record in records.values():
        if (mode=="exclude"):
            if (len(str(record.seq)) > length and not any(term in record.description for term in header_terms)):
                new_records[record.name] = record
        elif (mode=="include"):
            # print (record.description)
            # print (header_terms)
            # print (any(term in record.description for term in header_terms))
            # print (len(str(record.seq)))
            if (len(str(record.seq)) > length and any(term in record.description for term in header_terms)):
                new_records[record.name] = record

    return new_records

def addSpeciesFromHeader(records):
    """ Add species names that exist within the header of the seq ID i.e. like [Amazona aestiva] to a dictionary of
    seqRecords
    Add
    :param records: The dictionary to update with species names
    :return:
    """

    for record in records:
        if "[" in records[record].description:
            species_name = records[record].description.split(">")[0].split("[")[1].split("]")[0]
            records[record].annotations["Species"] = species_name
    return records

def getRecordsWithoutSpecies(records):
    """
    Get the entries in a seq records file that don't have a species annotated with them
    :param records: The records to check
    :return: The list of seq IDs missing species information
    """
    missing = []

    for record in records.values():
        if not "Species" in record.annotations:
            missing.append(record.id)

    return missing

def reportRecordsWithoutSpecies(records):
    """
    Print out statements about whether entries in a seq record file have annotated species
    :param records:
    :return:
    """
    missing = getRecordsWithoutSpecies(records)
    if missing:
        print("The following records do not have species annotated", [x for x in missing])
    else:
        print("All of the records have a species annotated")

def checkCDHitOutput(records, cdhit_output):
    """
    Reports on the species information that occurs within clusters in CD-hit output
    :param records:
    :param cdhit_output:
    :return:
    """

    if getRecordsWithoutSpecies(records):
        records = addSpeciesFromHeader(records)
        if getRecordsWithoutSpecies(records):
            print ("Some of your entries are missing their species annotation")
            return

    clusters = defaultdict(list)
    current = None

    # Read in and organise the CD-HIT output
    for record in records.values():

        for line in cdhit_output:
            if ("Cluster") in line:
                current = line

            if ("Cluster") not in line:

                clusters[current].append(re.sub(r"\s+", " ", line))

    # Only interested in clusters with more than one sequence
    for cluster in clusters.values():
        if len(cluster) > 1:
            for seq in cluster:
                trim = re.search('>(.*)\.\.\.', seq)
                seqID = trim.group(1).strip()
                if "*" in seq:
                    print ("Retained ", seqID, records[seqID].annotations["Species"])
                else:
                    print ("Removed", seqID, records[seqID].annotations["Species"])
            print ("")





def subset_records_with_regex(*header_terms, records, length=0, mode="exclude"):
    searchTerm = ""
    for term in header_terms:
        searchTerm += "\s" + term + "\d*-like|" + term + "\d*|"

    regex = re.compile(searchTerm[0:len(searchTerm) - 1])

    new_records = {}
    for record in records.values():
        if (mode=="exclude"):
            if (len(str(record.seq)) > length and not regex.search(record.description)):
                new_records[record.name] = record
        elif (mode=="include"):
            if (len(str(record.seq)) > length and regex.search(record.description)):
                new_records[record.name] = record

    return new_records

def exclude_character(records, character, mode="exclude"):
    new_records = {}
    for record in records.values():
        if (mode=="exclude"):
            if character not in record.seq:
                new_records[record.name] = record
        elif (mode == "incldue" ):
            if character in record.seq:
                new_records[record.name] = record
    return new_records



def build_species_count(*header_terms, records, length=0, mode="exclude"):
    species_count = defaultdict(list)
    for record in records.values():
        if (mode=="exclude"):
            if (len(str(record.seq)) > length and not any(term in record.description for term in header_terms)):
                if "[" in record.description:
                    name = record.description.split("[")[1].split("]")[0]
            else:
                print (record.description)
            species_count[name].append(record.id)
    return species_count




def count_ids(records, min=0):
    count = 0
    for record in records.values():
        if len(record) > min:
            count += len(record)
    return count

def get_species_names(records, min=0, counts=False):
    species_names = []
    for k,v in records.items():
        if len(v) > min:
            species_names.append(k) if not counts else species_names.append(k + " " + str(len(v)))
    return species_names

def plot_record_number(records, plotType, min=0):
    plot_records = {}
    # If we need to restrict the records to plot based on a minimum number
    if min > 0:

        for k,v in records.items():
            if len(v) >= min:
                plot_records[k] = v;
    else:
        plot_records = records

    if (plotType == "Bar"):

        data = [Bar(x=list(plot_records.keys()),
                    y= [len (x) for x in plot_records.values()])]

    elif (plotType == "Pie"):
        trace = go.Pie(labels=plot_records.keys(), values=plot_records.values())

        py.iplot([trace], filename='basic_pie_chart')



    else:
        print ("Invalid plot type selected. Choose 'Bar' or 'Pie'")

    return data

def map_dict_to_records(records, full_dict={}, unique=False):
    # If full_dict hasn't been provided, take records as the full dictionary
    if not full_dict:
        print ('empty dict')
        full_dict = records
    out_records = []
    for record in records.values():
        out_records.append(full_dict[record.id])
    return out_records

def map_species_dict_to_records(records, full_dict, unique=False):
    out_records = []
    for idlist in records.values():
        if unique:
            out_records.append(full_dict[idlist[0]])
        else:
            for id in idlist:
                out_records.append(full_dict[id])
    return out_records


def map_list_to_records(records, full_dict, unique=False):
    out_records = []
    for id in records:
        if unique and id in full_dict:
            out_records.append(full_dict[id])
        elif id in full_dict:
            out_records.append(full_dict[id])
    return out_records


def write_fasta(records, filename):
    with open(filename, "w") as handle:
            SeqIO.write(records, handle, "fasta")

def compare_fasta(records1, records2):
    for x in records1.keys():
        if x not in records2.keys():
            print (x)

def replaceWords(changes, file):
    for k,v in changes:
        file.replace(k,v)
    return file

def getDifferentIDs(records1, records2):
    ids1 = set(records1.keys())
    ids2 = set(records2.keys())
    return (ids1 - ids2)

def getDifferentRecords(records1, records2):
    ids = getDifferentIDs(records1, records2)
    differentDict = {k: records1[k] for k in ids}
    return differentDict

def subsetOnMotif(records, motif, retain_motif_seqs = True):
    """
    Subset a set of records based on the presence of a motif.

    :param motif: The motif to check for
    :param exclude: Whether we should retain sequences with the motif or exclude them
    :return: The subseted records based on the presence of the motif
    """
    subset_records = {}

    pattern=re.compile(motif)
    count = 1

    for record in records.values():
        if (len(pattern.findall(str(record.seq))) > 0 and retain_motif_seqs) or (len(pattern.findall(str(record.seq))) == 0 and not retain_motif_seqs):
            subset_records[record.id] = record
        else:
            print (count)
            print(record.id, pattern.findall(str(record.seq)))
            count +=1

    return subset_records



# full_record = SeqIO.to_dict(SeqIO.parse("files/candidates/regextest.fasta", "fasta"))
# full_record = SeqIO.to_dict(SeqIO.parse("/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and bacteria/180312_fifty_percent_identity/Alignment_exclusion/2U1_50_percent_tagged.fasta", "fasta"))
#
# motifSeqs = subsetOnMotif(full_record, '[FW][SGNH].[GD][^F][RKHPT][^P]C[LIVMFAP][GAD]', retain_motif_seqs=True)

# for record in motifSeqs.values():
#     print (record.name)
#
# write_fasta(motifSeqs.values(), "/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and bacteria/180312_fifty_percent_identity/Alignment_exclusion/2U1_50_percent_motif.fasta" )

# seq1 = "NPS"
# seq2 = "NASNAS"
# seq3 = "NAT"
#
# pattern = re.compile('N[^P][ST]')
#
# if len(pattern.findall(seq2)) > 0:
#     print ('found some')


# only_2U1_records = subset_records("2U1", "2U1-like", records=full_record, length=400, mode="include")
#
# print (only_2U1_records)

# test_records = subset_records_with_regex("2D", "2B", records=full_record, mode="include")
# for item in full_record:
#     print(full_record[item].description)
#
# print("---------------------------------")
#
# for item in test_records:
#     print(test_records[item].description)
