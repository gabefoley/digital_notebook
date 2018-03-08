from ete3 import NCBITaxa


Ltax = [561863, 333367, 518636, 1262999, 657322, 550540, 44012, 748224, 518636, 1309411]
ncbi = NCBITaxa()
tree = ncbi.get_topology(Ltax, intermediate_nodes=False)
print (tree)
# ncbi.annotate_tree(t, taxid_attr='name')

def get_species_to_rank_dict(tree, rank):
    rank_dict = {}
    for node in tree.iter_descendants("postorder"):
        if node.is_leaf():
            print (node)
            print (node.name)
            print(node.named_lineage)
            print(node.lineage)
            print (ncbi.get_rank(node.lineage))

            for key, val in ncbi.get_rank(node.lineage).items():
                if val == rank:
                    print ('got a rank')
                    rank_dict[node.name] = key
                    break
                else:
                    rank_dict[node.name] = ""


            print(node.sci_name)
            print(node.rank)
    return rank_dict
rank_dict = get_species_to_rank_dict(tree, "family")

print (rank_dict)
