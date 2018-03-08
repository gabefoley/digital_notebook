import ete3
import utilities

from ete3 import NCBITaxa

def updateTaxonomyDatabase():
    ncbi = NCBITaxa()
    ncbi.update_taxonomy_database()







# Pull the taxonomy information appended to the end of each ID (thanks SeqScrub!)

ncbi = NCBITaxa()

def get_species_to_taxonomy_dict(tree):
    taxonomySet = set()
    for leaf in tree.get_leaf_names():
        species = leaf.split("_")[-1]
        if species[0].isupper():
            taxonomySet.add(species)

    # Map the species names to taxon IDs
    taxonomy_dict = ncbi.get_name_translator(taxonomySet)

    for key, val in taxonomy_dict.items():
        taxonomy_dict[key] = str(val[0])

    #TODO: Fix this so that the dict values aren't in a list
    return taxonomy_dict

def get_id_to_taxonomy_dict(tree, seq_type="protein"):
    seq_ids = []
    for leaf in tree.get_leaf_names():
        seq_id = leaf.split(" ")[0]
        seq_ids.append(seq_id)

    # Map the seq IDs to taxon IDs
    taxonomy_dict = utilities.buildTaxonomyDict(seq_ids, seq_type)
    return taxonomy_dict


# print (taxonomySet)
# name2taxid = ncbi.get_name_translator(list(taxonomySet))
# print (name2taxid)

def create_taxon_tree_from_species(tree, taxonomy_dict):

    taxon_tree = tree.copy()

    for leaf in taxon_tree:
        seq_id = leaf.name.split("_")[-1]
        leaf.name = taxonomy_dict[seq_id]

    return taxon_tree

def create_taxon_tree_from_id(tree, taxonomy_dict):

    taxon_tree = tree.copy()

    for leaf in taxon_tree:
        seq_id = leaf.name.split(" ")[0]
        leaf.name = taxonomy_dict[seq_id]


    return taxon_tree

def annotate_nodes_based_on_taxonomy(tree):
    """
    Takes a tree with leaves already labelled with taxonomic ID and labels internal nodes with the union of child names
    :param taxon_tree: The tree labelled with taxonomic ID
    :return: A tree with interally labelled nodes based on union of taxonomic ID
    """

    taxon_tree = tree.copy()

    for node in taxon_tree.traverse(strategy="postorder"):

        children_list = []

        if not node.is_leaf():
            for desc in node.children:

                for child in str(desc.name).split(","):
                    children_list.append(child)

            node.name = ",".join(str(x) for x in set(children_list))

    return taxon_tree

def get_ncbi_tree(taxonomy_dict):
    """
    Create a tree based on the NCBI taxonomy database from a set of taxonomic IDs
    :param taxonomy_dict: The taxonomy dictionary that contains the taxonomic IDs as its values
    :return: A tree based on the NCBI taxonomy database
    """
    ncbi_tree = ncbi.get_topology([x for x in taxonomy_dict.values()])
    return ncbi_tree


# print (latestTree.get_ascii(show_internal=True))
# print (node2labels)
# print (latestTree.write(is_leaf_fn=collapsed_leaf))


# print (ncbi_tree.get_ascii(attributes=["sci_name", "rank"]))
# print (collapsed_tree)
# print (latestTree)

def collapse_tree(tree):


    def collapsed_leaf(node):
        if len(node2labels[node]) == 1:
            return True
        else:
            return False

    node2labels = tree.get_cached_content(store_attr="name")
    collapsed_tree = ete3.Tree(tree.write(is_leaf_fn=collapsed_leaf))



    return collapsed_tree

def get_duplicate_leaves(tree):
    """
    Takes a phylogenetic tree and returns a list of duplicated leaf names (if any)
    :param tree: The phylogenetic tree to check
    :return: A list of duplicated leaf names if there are duplicates, otherwise an empty list
    """

    leaf_names = []
    duplicate_leaf_names = set()
    for leaf in tree.get_leaves():
        if leaf.name in leaf_names:
            duplicate_leaf_names.add(leaf.name)
        else:
            leaf_names.append(leaf.name)

    return list(duplicate_leaf_names)

def get_species_to_rank_dict(tree, rank):
    rank_dict = {}
    for node in tree.iter_descendants("postorder"):
        if node.is_leaf():

            for key, val in ncbi.get_rank(node.lineage).items():
                if val == rank:
                    rank_dict[node.name] = str(key)
                    break
                else:
                    rank_dict[node.name] = "XXXX"
    return rank_dict

def compare_tree_with_taxon_tree(*args, filepath, id="True", filepath_to_save_dict="", filepath_to_load_dict=""):

    tree = ete3.Tree(filepath, format=1)

    # If user wants to load the taxonomy dictionary from file
    if filepath_to_load_dict:
        taxonomy_dict = utilities.open_python_object(filepath_to_load_dict)

    if id:
        # If we don't already have a taxonomy dict, create one
        if 'taxonomy_dict' not in locals():
            taxonomy_dict = get_id_to_taxonomy_dict(tree)
        taxon_tree = create_taxon_tree_from_id(tree, taxonomy_dict)

    else:
        # If we don't already have a taxonomy dict, create one
        if  'taxonomy_dict' not in locals():
            taxonomy_dict = get_species_to_taxonomy_dict(tree)
        taxon_tree = create_taxon_tree_from_species(tree, taxonomy_dict)

    # If user wants to save the taxonomy dictionary to file
    if filepath_to_save_dict:
        utilities.save_python_object(taxonomy_dict, filepath_to_save_dict)


    # taxon_tree = annotate_nodes_based_on_taxonomy(taxon_tree)

    # collapsed_tree = collapse_tree(taxon_tree)
    print ("\nThis is the original tree")
    print (taxon_tree)

    # print ("The following taxon ids appear in multiple clades of the tree", duplicate_leaves)
    ncbi_tree = get_ncbi_tree(taxonomy_dict)


    print ("\nThis is the NCBI species tree\n")
    print (ncbi_tree)

    for rank in args:

        rank_dict = get_species_to_rank_dict(ncbi_tree, rank)
        ncbi_rank_tree = create_taxon_tree_from_id(ncbi_tree, rank_dict)

        rank_tree = create_taxon_tree_from_id(taxon_tree, rank_dict)




        print("\nThis is the original tree based on %s \n" % rank)
        print(rank_tree)

        print ("\nThis is the NCBI tree based on %s \n" % rank)
        print(ncbi_rank_tree)

        annotated_rank_tree = annotate_nodes_based_on_taxonomy(rank_tree)
        collapsed_rank_tree = collapse_tree(annotated_rank_tree)


        print ("\nThis is the collapsed original tree based on %s \n" % rank)
        print (collapsed_rank_tree)

        annotated_ncbi_rank_tree = annotate_nodes_based_on_taxonomy(ncbi_rank_tree)
        collapsed_ncbi_rank_tree = collapse_tree(annotated_ncbi_rank_tree)

        print ("\nThis is the collapsed NCBI tree based on %s \n" % rank)
        print (collapsed_ncbi_rank_tree)


        # Check if the collapsed rank tree has duplicate leaves (we won't be able to perform Robinson Foulds if it does
        duplicate_leaves = get_duplicate_leaves(collapsed_rank_tree)

        if not duplicate_leaves:

            rf, max_rf, common_leaves, parts_t1, parts_t2, set1, set2 = collapsed_rank_tree.robinson_foulds(collapsed_ncbi_rank_tree, unrooted_trees=True)
            # print ("\nRF distance is %s over a total of %s" %(rf, max_rf))
            # print ("Partitions in tree1:", parts_t1)
            # print ("Partitions in tree2:", parts_t2)
            print ("\nPartitions in tree2 that were not found in tree1:", parts_t1 - parts_t2)
            print ("Partitions in tree1 that were not found in tree2:", parts_t2 - parts_t1)

        else:
            print ("\nThese leaves were duplicated in the original collapsed %s tree" % (rank))
            print (duplicate_leaves)



# print (latestTree)
# ncbi_tree.write(format=3, outfile="species_tree.nwk")
# latestTree.write(format=1, outfile="latest_tree.nwk")
# latestTree = ete3.Tree("files/2U1_tree_rerooted_removed.nwk", format=1)
# latestTree = ete3.Tree("files/2U1_removed.nwk", format=1)

# test = ["files/taxon_annotate"]
test = ["files/test"]
filepaths = ["files/2U1_complete_sstrimmed_cleaned_phyml", "files/2U1_N_terminal_intact_trimmed_cleaned_phyml", "files/2U1_N_terminal_intact_insertions_removed_trimmed_cleaned_phyml" ]
for path in test:
    tax_dict_path = path + "_taxon_dict.obj"
    compare_tree_with_taxon_tree("family", filepath=path + ".nwk", id=True, filepath_to_load_dict=tax_dict_path)