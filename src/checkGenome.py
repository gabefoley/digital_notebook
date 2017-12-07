from Bio import SeqIO, Entrez
from lxml import etree as ET

from Bio.SeqFeature import SeqFeature, FeatureLocation
my_seq_feature = SeqFeature(FeatureLocation(50,100),strand=+1)

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib.units import cm

Entrez.email = "gabriel.foley@uqconnect.edu.au"


def check_genomic_location(records, min=0, visualise=None, file_path=None):
    for species, ids in records.items():
        if len(ids) > min:
            records = ""
            recordDict = {}

            for record in ids:
                records += record + " OR "

            # print (species)
            # records = records[0:-4]

                search = Entrez.esearch(term=record, db="gene", retmode="fasta", rettype="acc")
                result = Entrez.read(search)
                # print ('progress')
                # print (result)
                if len(result['IdList']) > 0:
                    gene_ids = result['IdList'][0]
                    # print ('here are gene ids')
                    # print (gene_ids)
                    recordDict[record] = gene_ids

            print (recordDict)
            # print (len(gene_ids))

            if len(recordDict) > 1:

                if not visualise:

                    for protein_id, gene_id in recordDict.items() :
                        handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
                        record = handle.read()
                        xml_parsed = ET.fromstring(record)

                        start = xml_parsed.xpath("/Entrezgene-Set/Entrezgene/Entrezgene_track-info/Gene-track/Gene-track_geneid[contains(., '" + gene_id + "')]/../../../Entrezgene_locus/Gene-commentary/Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval/Seq-interval_from/text()")
                        finish = xml_parsed.xpath("/Entrezgene-Set/Entrezgene/Entrezgene_track-info/Gene-track/Gene-track_geneid[contains(., '" + gene_id + "')]/../../../Entrezgene_locus/Gene-commentary/Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval/Seq-interval_to/text()")
                        chromosome = xml_parsed.xpath("/Entrezgene-Set/Entrezgene/Entrezgene_track-info/Gene-track/Gene-track_geneid[contains(., '" + gene_id + "')]/../../../Entrezgene_source/BioSource/BioSource_subtype/SubSource/SubSource_name/text()")

                        if file_path == None:
                            print()
                            print ("Protein id is %s" % (protein_id))
                            print ("Gene"
                                   " id is %s " % (gene_id))
                            print("Gene region starts at %s" % (start[0]))
                            print ("Gene region ends at %s" % (finish[0]))
                            print ("Gene region is %s nucleotides long " % (int(finish[0]) - int(start[0])))
                            if (chromosome[0] == "Un"):
                                print ("Chromosome is unassigned")
                            else:
                                print ("On chromosome %s" % (chromosome[0]))

                            print ("-------------------------------------------------------------------")

                        else:
                            with open(file_path + " " + species, "a+") as text_file:
                                text_file.write("Protein id is %s \n" % (protein_id))
                                text_file.write("Gene"
                                                " id is %s \n" % (gene_id))
                                text_file.write("Gene wregion starts at %s \n" % (start[0]))
                                text_file.write("Gene region ends at %s \n" % (finish[0]))
                                text_file.write("Gene region is %s nucleotides long \n" % (int(finish[0]) - int(start[0])))
                                if (chromosome[0] == "Un"):
                                    text_file.write("Chromosome is unassigned \n")
                                else:
                                    text_file.write("On chromosome %s \n" % (chromosome[0]))

                                text_file.write("-------------------------------------------------------------------\n")

                elif (visualise == "circular" or "linear"):
                    drawGenome(species, recordDict, visualise)

                else:
                    print ("-visualise flag should be either 'circular' or 'linear'")

def drawGenome(species, recordDict, visualise):

    locations = []

    gdd = GenomeDiagram.Diagram('Test Diagram')
    gdt_features = gdd.new_track(1, greytrack=False)
    gds_features = gdt_features.new_set()

    for protein_id, gene_id in recordDict.items():
        handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
        record = handle.read()
        xml_parsed = ET.fromstring(record)

        start = xml_parsed.xpath(
            "/Entrezgene-Set/Entrezgene/Entrezgene_track-info/Gene-track/Gene-track_geneid[contains(., '" + gene_id + "')]/../../../Entrezgene_locus/Gene-commentary/Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval/Seq-interval_from/text()")
        finish = xml_parsed.xpath(
            "/Entrezgene-Set/Entrezgene/Entrezgene_track-info/Gene-track/Gene-track_geneid[contains(., '" + gene_id + "')]/../../../Entrezgene_locus/Gene-commentary/Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval/Seq-interval_to/text()")
        if (start and finish):

            locations.append(int(start[0]))
            locations.append(int(finish[0]))
            feature = SeqFeature(FeatureLocation(int(start[0]), int(finish[0])), strand=+1)
            gds_features.add_feature(feature, name=protein_id + "(" + gene_id + ")", label=True)

    gdd.draw(format=visualise, pagesize=(15*cm,20*cm), fragments=1,
                 start=0, end=max(locations))
    gdd.write(species + ".pdf", "pdf")