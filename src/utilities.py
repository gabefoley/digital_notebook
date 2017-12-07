from Bio import SeqIO, SeqRecord

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

# loadSequences("files/candidates/Amazona_aestiva.fasta")

# loadSequences( "files/candidates/Andrias_davidianus.fasta",
#             "files/candidates/Marmota_marmota.fasta", "files/candidates/Poecilia_reticulata.fasta",
#             "files/candidates/Python_bivittatus.fasta")

def saveHeaderTerms(header_terms, filepath):
    header_string = " ".join(header_terms)
    print (header_string)
    print (type(header_string))
    with open(filepath, "w") as text_file:
        text_file.write(header_string)

def loadHeaderTerms(filepath):

    header_terms = []
    with open(filepath, "r") as text_file:
        for item in text_file.split():
            header_terms.append(item)
    return header_terms


def addHeaderTerms(sender, header_terms):
    for item in sender.value.split():
        header_terms.append(item)
