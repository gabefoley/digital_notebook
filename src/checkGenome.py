from Bio import SeqIO, Entrez
from lxml import etree as ET

from Bio.SeqFeature import SeqFeature, FeatureLocation
my_seq_feature = SeqFeature(FeatureLocation(50,100),strand=+1)

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib.units import cm

Entrez.email = "gabriel.foley@uqconnect.edu.au"


def check_genomic_location(records, min=0, visualise=None):
    for species, ids in records.items():
        if len(ids) > min:
            records = ""

            for record in ids:
                records += record + " OR "

            print (species)
            records = records[0:-4]

            search = Entrez.esearch(term=records, db="gene", retmode="fasta")
            result = Entrez.read(search)
            gene_ids = result['IdList']
            print (gene_ids)
            print (len(gene_ids))

            if len(gene_ids) > 1:

                if not visualise:

                    for gene_id in gene_ids:
                        handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
                        record = handle.read()
                        xml_parsed = ET.fromstring(record)

                        start = xml_parsed.xpath("/Entrezgene-Set/Entrezgene/Entrezgene_track-info/Gene-track/Gene-track_geneid[contains(., '" + gene_id + "')]/../../../Entrezgene_locus/Gene-commentary/Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval/Seq-interval_from/text()")
                        finish = xml_parsed.xpath("/Entrezgene-Set/Entrezgene/Entrezgene_track-info/Gene-track/Gene-track_geneid[contains(., '" + gene_id + "')]/../../../Entrezgene_locus/Gene-commentary/Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval/Seq-interval_to/text()")
                        chromosome = xml_parsed.xpath("/Entrezgene-Set/Entrezgene/Entrezgene_track-info/Gene-track/Gene-track_geneid[contains(., '" + gene_id + "')]/../../../Entrezgene_source/BioSource/BioSource_subtype/SubSource/SubSource_name/text()")


                        print()
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

                elif (visualise == "circular" or "linear"):
                    drawGenome(species, gene_ids, visualise)

                else:
                    print ("-visualise flag should be either 'circular' or 'linear'")

def drawGenome(species, gene_ids, visualise):

    locations = []

    gdd = GenomeDiagram.Diagram('Test Diagram')
    gdt_features = gdd.new_track(1, greytrack=False)
    gds_features = gdt_features.new_set()

    for gene_id in gene_ids:
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
            gds_features.add_feature(feature, name=gene_id, label=True)

    gdd.draw(format=visualise, pagesize=(15*cm,4*cm), fragments=1,
                 start=0, end=max(locations))
    gdd.write(species + ".pdf", "pdf")