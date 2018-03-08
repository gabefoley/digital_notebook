from Bio import SeqIO, SeqRecord
from Bio import Entrez
import pandas as pd
import random
import os
import pickle

Entrez.email = "gabriel.foley@uqconnect.edu.au"


def loadSequences(*args, unique=True):
    print ("lets load")
    '''
    Join multiple sequence files together
    :param args: The sequence files to open
    :param unique: Boolean representing whether we only want unique sequences
    :return:
    '''
    for file in args:
        print (file)
        current_handle2 = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
        print (len(current_handle2))


def saveIDs(*args, percent_identity="", output_dir="", concatenate=False):
    '''
    Load in a list of IDs and retreive the sequences for them
    :param args: The lists to open
    :param filepath: Optional filepath to write file to
    :return: SeqRecord object if no file path specified, otherwise nothing
    '''

    id_list = []
    output_id = randstring(5)


    for file in args:
        outpath = file[file.rindex('/')+1:].replace(".", "_output.")



        # df = pd.read_csv(file, delimiter='\t', header=None)
        df = pd.read_csv(file, delimiter='\t', comment='#', names= ["query acc.ver", "subject acc.ver", "% identity",
                                                                    "alignment length", "mismatches", "gap opens",
                                                                    "q. start", "q. end", "s. start", "s. end",
                                                                    "evalue", "bit score"])

        #If user has set a minimum percent identity then filter by that
        if percent_identity:
            id_list += df.loc[df['% identity'] >= percent_identity, 'subject acc.ver'].drop_duplicates().tolist()

        # No minimum percent identity so take all records
        else:
            id_list = df['subject acc.ver'].drop_duplicates().tolist()

        # If we're not concatenating lets write this out to the
        if not concatenate:
            # Remove any existing previous file at this location
            removeFile(output_dir + outpath)

            # Write out the sequences to the output path
            getIDs(id_list, output_dir + outpath)
            id_list = []

    if concatenate:
        # Get rid of any duplicates in id_list (i.e. hits that appeared in more than one input file)
        id_list = list(set(id_list))
        getIDs(id_list, output_dir + output_id +"_concatenated_output.fasta")

def loadIDs(*args, percent_identity=""):
    id_list = []

    for file in args:

        # df = pd.read_csv(file, delimiter='\t', header=None)
        df = pd.read_csv(file, delimiter='\t', comment='#', names= ["query acc.ver", "subject acc.ver", "% identity",
                                                                    "alignment length", "mismatches", "gap opens",
                                                                    "q. start", "q. end", "s. start", "s. end",
                                                                    "evalue", "bit score"])

        #If user has set a minimum percent identity then filter by that
        if percent_identity:
            id_list += df.loc[df['% identity'] >= percent_identity, 'subject acc.ver'].drop_duplicates().tolist()

        # No minimum percent identity so take all records
        else:
            id_list = df['subject acc.ver'].drop_duplicates().tolist()

    #Remove any duplicates in id_list (i.e. hits that appeared in more than one input file)
    id_list = list(set(id_list))
    handle = getIDs(id_list)
    return handle



def getIDs(id_list, filepath=""):
    for i in range(0, len(id_list), 500):

        final = i + 500 if (i + 500 < len(id_list)) else len(id_list) + 1

        handle = Entrez.efetch(db="protein", rettype="fasta", retmode="text", \
                               id=id_list[i:final])

        if filepath:
            with open(filepath, "a") as outfile:
                for seq_record in SeqIO.parse(handle, "fasta"):
                    outfile.write(">" + seq_record.description +"\n")
                    outfile.write(str(seq_record.seq + "\n"))
                handle.close()

        else:
            return SeqIO.parse(handle, "fasta")


def saveHeaderTerms(header_terms, filepath):
    header_string = " ".join(header_terms)
    print (header_string)
    print (type(header_string))
    with open(filepath, "w") as text_file:
        text_file.write(header_string)

def loadHeaderTerms(filepath):

    header_terms = []
    with open(filepath, "r") as text_file:
        for item in text_file.read().split():
            header_terms.append(item)
    return header_terms


def addHeaderTerms(sender, header_terms):
    for item in sender.value.split():
        header_terms.append(item)


def buildTaxonomyDict(seq_ids, seq_type="protein"):
    """
    Take a list of sequence ids and return a dictionary mapping those ids to their taxonomic ID

    :param seq_ids: List of sequence ids
    :param seq_type: Type of sequence / which database to query
    :return: Dictionary mapping seq ID to taxonomic ID
    """

    taxonomy_dict = {}

    for seq_id in seq_ids:
        handle = Entrez.elink(dbfrom=seq_type, db="taxonomy", id=seq_id)
        records = Entrez.read(handle)
        print (records)
        taxonomy_dict[seq_id] = records[0]["LinkSetDb"][0]["Link"][0]['Id']

    return taxonomy_dict

def randstring(length=10):
    valid_letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    return ''.join((random.choice(valid_letters) for i in range(length)))

def removeFile(*args):
    """
    Remove files in the list from the directory

    :param args: Files to remove
    :return:
    """
    for arg in args:
        print (arg)
        if (os.path.exists(arg)):
            print ("exists")
            os.remove(arg)

def save_python_object(object, filepath):
    file_path = open(filepath, 'wb')
    pickle.dump(object, file_path)

def open_python_object(filepath):
        file = open(filepath, 'rb')
        object = pickle.load(file)
        return object

# loadIDs("/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/id_list.txt", percent_identity=30, filepath="/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/BLAST_unique.fasta" )
# saveIDs("/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants and fungi/Andrias_davidianus.txt",
#         "/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants and fungi/Python_bivittatus.txt",
#         percent_identity=83, output_dir = "/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants and fungi/",  concatenate=True )

directory = "/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and bacteria/"
saveIDs(directory + "Amazona_aestiva.txt", directory + "Andrias_davidianus.txt", directory + "Marmota_marmota.txt",  directory + "Poecilia_reticulata.txt", directory + "Python_bivittatus.txt", percent_identity=50, output_dir = directory,  concatenate=True)
