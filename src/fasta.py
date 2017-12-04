from Bio import SeqIO, Entrez
from collections import defaultdict
from lxml import etree as ET
import plotly
plotly.tools.set_credentials_file(username='gabefoley', api_key='xS8qT0kIbIKDWt0BalOd')
import plotly.plotly as py
from plotly.graph_objs import *


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

def build_species_count(*exclude_on, records, length=0 ):
    species_count = defaultdict(list)
    for record in records.values():
        if (len(str(record.seq)) > length and not any(exclude in record.description for exclude in exclude_on)):
            if "[" in record.description:
                name = record.description.split("[")[1].split("]")[0]
            else:
                print (record.description)
            # print (record.annotations["organism"])
            species_count[name].append(record.id)
    return species_count




def count_ids(records, min=0):
    count = 0
    for record in records.values():
        if len(record) > min:
            count += len(record)
    return count

def get_species_names(records, min=0, counts=False):
    for k,v in records.items():
        if len(v) > min:
            print (k) if not counts else print (k, len(v))

def plot_record_number(records, min=0):
    plot_records = {}
    # If we need to restrict the records to plot based on a minimum number
    if min > 0:

        for k,v in records.items():
            if len(v) >= min:
                plot_records[k] = v;
    else:
        plot_records = records

    data = [Bar(x=list(plot_records.keys()),
                    y= [len (x) for x in plot_records.values()])]

    return data


def map_ids_to_records(records, full_dict, unique=False):
    out_records = []
    for idlist in records.values():
        if unique:
            out_records.append(full_dict[idlist[0]])
        else:
            for id in idlist:
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
