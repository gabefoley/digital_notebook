from Bio import SeqIO, Entrez

Entrez.email = "gabriel.foley@uqconnect.edu.au"
homo_sapiens = SeqIO.to_dict(SeqIO.parse("files/exons/homo_sapiens.fasta", "fasta"))
strongylocentrotus_purpuratus = SeqIO.to_dict(SeqIO.parse("files/exons/strongylocentrotus_purpuratus.fasta", "fasta"))
all = SeqIO.to_dict(SeqIO.parse("files/exons/all.fasta", "fasta"))

def mapExons (records):
    for record in records:

        handle = Entrez.elink(dbfrom="protein", db='gene', id=records[record].id)
        #
        mapping = Entrez.read(handle)

        for term in mapping:
            geneID = (term['LinkSetDb'][0]['Link'][0]['Id'])

        handle = Entrez.efetch(db="protein", id=records[record].id, rettype="gb", retmode="xml")
        full_record = Entrez.read(handle)

        for term in full_record:
            # print (term.keys
            # print (term)
            # print (term['GBSeq_feature-table'])
            # #
            # print (term['GBSeq_feature-table'][3])
            refseq = (term['GBSeq_source-db'].split("accession")[1].strip())
            # print (term.keys())
            # for item in term:
            #     print (item)
            #     print ("******8")
            #     print (term[item])
        print ("The original protein record is %s which has a gene ID %s and a RefSeq ID %s " % (record, geneID, refseq))

        handle = Entrez.efetch(db="gene", id=geneID, rettype="gb", retmode="xml")

        gene_record = Entrez.read(handle)

        # print (gene_record)
        for term in gene_record:

            for j in term['Entrezgene_comments']:
                if 'Gene-commentary_comment' in j:
                    print (j)
                    print (j.keys())
                    print(j['Gene-commentary_comment'][0])
                    # print (j['Gene-commentary_products'][0]['Gene-commentary_accession'])
                    # [0]['Gene-commentary_products'][0]['Gene-commentary_accession'])
            # print (term)
            # print (term['Entrezgene_comments'])
            # # print (term['Entrezgene_comments']['Gene-commentary_comment'])

            # print (term['Entrezgene_comments'][2]['Gene-commentary_comment'][0]['Gene-commentary_products'][0]['Gene-commentary_accession'])
            print (term['Entrezgene_comments'][2]['Gene-commentary_comment'][0]['Gene-commentary_products'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_from'])
            print (term['Entrezgene_comments'][2]['Gene-commentary_comment'][0]['Gene-commentary_products'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_to'])



mapExons(strongylocentrotus_purpuratus)
mapExons(all)