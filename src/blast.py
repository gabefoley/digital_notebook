from Bio.Blast.Applications import NcbipsiblastCommandline, NcbitblastnCommandline
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW

import sys
import subprocess


def psi_blast(db_file, query_file, eval_num, iter_num):

    # fasta_string = open(queryFile)
    # result_handle = NCBIWWW.qblast("blastp", "nr", fasta_string)
    # print (result_handle)
    # for record in result_handle:
    #     print (record)

    # queryFile = open(query).read()

    psi_cline = NcbipsiblastCommandline('psiblast', db=db_file, query=query_file, evalue=eval_num,
                                        out="files/psi3.xml", outfmt=5, out_pssm="_pssm", remote="yes")
    p = subprocess.Popen(str(psi_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         shell=(sys.platform != "win32"))
    print(psi_cline)
    print(p)


def get_blast_info(xml_file):
    """
    Get the location and sequence of the top scoring hits from a BLAST results
    :param xml_file: BLAST records
    :return:
    """
    blast_info = {}
    blast_parser = NCBIXML.parse(xml_file)
    for record in blast_parser:
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                blast_info["location"] = str(hsp.sbjct_start) + ":" + str(hsp.sbjct_end)
                blast_info["sequence"] = hsp.sbjct
                break
    return blast_info


# psi_blast("nr", "files/camelus_dromedarius.fasta", 0.00005, 5)
# get_blast_info(open("files/psi.xml"))
