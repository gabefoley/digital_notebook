from Bio import SeqIO, SeqRecord, Entrez
import re

Entrez.email = "gabefoley@gmail.com"
def getCommonNameFromSpecies(*species, seqtype="protein"):
    pass

def getCommonNameFromID(*ids, seq_type="protein"):
    taxonomy_dict = {}

    for seq_id in ids:
        handle = Entrez.efetch(db="protein", id=seq_id, retmode="xml")

        records = Entrez.read(handle)
        lineage = records[0]['GBSeq_taxonomy']
        name = records[0]['GBSeq_source']
        if ("(" in name):
            name = re.search(r'\((.*?)\)', name).group(1)
            print (lineage)
            print (name)
        else:
            print (lineage)
            # print (name + " ** NOT COMMON NAME **")




getCommonNameFromID('PIO41157.1', 'XP_013766542.1', 'PIK57532.1', 'ELW64418.1', 'KYO44822.1', 'XP_004941452.1', 'XP_006768464.1', 'XP_023361104.1', 'OCA35976.1', 'EDL30916.1', 'XP_016404253.1', 'XP_021107397.1', 'XP_014379039.1', 'XP_021107396.1', 'KTF79201.1', 'XP_016094679.1', 'XP_006123186.1', 'XP_014801266.1', 'XP_002605102.1', 'EPY88763.1', 'XP_006520454.1', 'PIK38219.1', 'XP_006825012.1', 'XP_013865517.1')