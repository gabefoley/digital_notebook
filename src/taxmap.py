

from Bio import SeqIO, Entrez
Entrez.email = "gabriel.foley@uqconnect.edu.au"




tax = buildTaxonomyDict(["XP_010584537.1"], "protein")


print (tax)