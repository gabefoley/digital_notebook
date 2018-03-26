import ete3
import utilities


def update_taxonomy_database():
    """
    Function to update the NCBI taxonomy database
    """
    ncbi = ete3.NCBITaxa()
    ncbi.update_taxonomy_database()

ncbi = ete3.NCBITaxa()


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
    """
    Create a dictionary that maps sequence ID to taxonomy ID
    :param tree: The tree to extract the sequence IDs from
    :param seq_type: The type of sequences
    :return: The dictionary mapping sequence ID to taxonomy ID
    """
    seq_ids = []
    for leaf in tree.get_leaf_names():
        seq_id = leaf.split("|")[0]
        seq_ids.append(seq_id)

    # Map the seq IDs to taxon IDs
    taxonomy_dict = utilities.build_taxonomy_dict(seq_ids, seq_type)
    return taxonomy_dict


def create_taxon_tree_from_species(tree, taxonomy_dict):

    taxon_tree = tree.copy()

    for leaf in taxon_tree:
        seq_id = leaf.name.split("_")[-1]
        leaf.name = taxonomy_dict[seq_id]

    return taxon_tree

def create_taxon_tree_from_id(tree, taxonomy_dict):

    taxon_tree = tree.copy()

    for leaf in taxon_tree:
        seq_id = leaf.name.split("|")[0]
        if seq_id in taxonomy_dict:
            leaf.name = taxonomy_dict[seq_id]
        elif seq_id != "XXXX" and seq_id[0:len(seq_id) - 3] != "XXXX":
            leaf.name = taxonomy_dict[seq_id[0:len(seq_id) - 3]] + utilities.randstring(3) # If we've assigned a random
            # string then remove it so that we can map back to the taxonomy dictionary

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

def change_duplicate_leaves(tree):
    """
    Takes a phylogenetic tree and changes duplicated leaf names (if any)
    :param tree: The phylogenetic tree to check
    :return: The changed tree
    """

    leaf_names = []
    for leaf in tree.get_leaves():
        if leaf.name in leaf_names:
            leaf.name = leaf.name + utilities.randstring(3)
        else:
            leaf_names.append(leaf.name)

    return tree

def get_unique_tree(tree, rank):
    """
    Check if a tree has duplicate leaves and change tree into a unique tree
    :param tree:
    :return:
    """


    duplicate_leaves = get_duplicate_leaves(tree)

    if duplicate_leaves:
        print("\nThese leaves were duplicated in the original collapsed %s tree" % (rank))
        print(duplicate_leaves)


    unique_tree = change_duplicate_leaves(tree)
    return unique_tree


def check_robinson_foulds(tree1, tree2, tree1_name="tree1", tree2_name="tree2"):
    rf, max_rf, common_leaves, parts_t1, parts_t2, set1, set2 = tree1.robinson_foulds(
        tree2, unrooted_trees=True)
    print ("\nRF distance between %s and %s is %s over a total of %s" %(tree1_name, tree2_name, rf, max_rf))
    print("Partitions in %s that were not found in %s: %s" % (tree1_name, tree2_name, parts_t2 - parts_t1))
    print("\nPartitions in %s that were not found in %s: %s" % (tree2_name, tree1_name, parts_t1 - parts_t2))


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

def get_rank_to_lineage_names(tree):
    """
    Create a dictionary that maps taxonomic rank to the actual lineage names
    :param tree:
    :return:
    """
    lineage_dict = {}
    for node in tree.iter_descendants("postorder"):
        if node.is_leaf():
            lineage_names = ""
            for taxid in node.lineage:
                if ncbi.get_rank(node.lineage)[taxid] != "no rank":

                    for key, val in ncbi.get_taxid_translator([taxid]).items():
                        if key != 1:  # Don't include the root
                            lineage_names += val + "_"
            lineage_names = lineage_names[0:len(lineage_names)-1] # Remove the final underscore
            lineage_dict[node.name] = lineage_names
    return lineage_dict

def get_rank_to_taxon_name(tree):
    """
    Create a dictionary that maps taxonomic rank to the actual taxon name
    :param tree:
    :return:
    """
    taxon_dict = {}
    for node in tree.iter_descendants("postorder"):
        if node.is_leaf():
            if node.name != "XXXX" and node.name[0:len(node.name)-3] != "XXXX":
                for taxon_name in ncbi.get_taxid_translator([int(node.name)]).values():
                    taxon_dict[node.name] = taxon_name
    return taxon_dict

def compare_tree_with_taxon_tree(*args, filepath, id="True", print_trees = False, filepath_to_save_dict="", filepath_to_load_dict="", outpath=""):

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
        if 'taxonomy_dict' not in locals():
            taxonomy_dict = get_species_to_taxonomy_dict(tree)
        taxon_tree = create_taxon_tree_from_species(tree, taxonomy_dict)

    # If user wants to save the taxonomy dictionary to file
    if filepath_to_save_dict:
        utilities.save_python_object(taxonomy_dict, filepath_to_save_dict)


    # taxon_tree = annotate_nodes_based_on_taxonomy(taxon_tree)

    # collapsed_tree = collapse_tree(taxon_tree)
    if print_trees:
        print ("\nThis is the original tree")
        print (taxon_tree)

    # print ("The following taxon ids appear in multiple clades of the tree", duplicate_leaves)
    ncbi_tree = get_ncbi_tree(taxonomy_dict)

    if print_trees:
        print ("\nThis is the NCBI species tree\n")
        print (ncbi_tree)

    check_robinson_foulds(taxon_tree, ncbi_tree, tree1_name="Original", tree2_name="NCBI")




    for rank in args:

        rank_dict = get_species_to_rank_dict(ncbi_tree, rank)
        ncbi_rank_tree = create_taxon_tree_from_id(ncbi_tree, rank_dict)

        rank_tree = create_taxon_tree_from_id(taxon_tree, rank_dict)

        if print_trees:
            print("\nThis is the original tree based on %s \n" % rank)
            print(rank_tree)

            print ("\nThis is the NCBI tree based on %s \n" % rank)
            print(ncbi_rank_tree)

        annotated_rank_tree = annotate_nodes_based_on_taxonomy(rank_tree)
        collapsed_rank_tree = collapse_tree(annotated_rank_tree)


        if print_trees:
            print ("\nThis is the collapsed original tree based on %s \n" % rank)
            print (collapsed_rank_tree)

        annotated_ncbi_rank_tree = annotate_nodes_based_on_taxonomy(ncbi_rank_tree)
        collapsed_ncbi_rank_tree = collapse_tree(annotated_ncbi_rank_tree)

        if print_trees:
            print ("\nThis is the collapsed NCBI tree based on %s \n" % rank)
            print (collapsed_ncbi_rank_tree)

        # Get a unique version of the collapsed rank tree (we won't be able to perform Robinson Foulds if it isn't)

        unique_collapsed_rank_tree = get_unique_tree(collapsed_rank_tree, rank)
        check_robinson_foulds(unique_collapsed_rank_tree, collapsed_ncbi_rank_tree, tree1_name="Supplied tree on basis of " + rank, tree2_name="NCBI tree on basis of " + rank)




    # Create lineage dictionaries
    lineage_dict = get_rank_to_lineage_names(ncbi_tree)
    lineage_tree = create_taxon_tree_from_id(taxon_tree, lineage_dict)
    ncbi_lineage_tree = create_taxon_tree_from_id(ncbi_tree, lineage_dict)

    lineage_rank_dict = get_rank_to_taxon_name(collapsed_ncbi_rank_tree)
    lineage_rank_tree = create_taxon_tree_from_id(collapsed_rank_tree, lineage_rank_dict)
    unique_lineage_rank_tree = create_taxon_tree_from_id(unique_collapsed_rank_tree, lineage_rank_dict)

    ncbi_lineage_rank_tree = create_taxon_tree_from_id(collapsed_ncbi_rank_tree, lineage_rank_dict)
    # unique_ncbi_lineage_rank_tree = create_taxon_tree_from_id(unique_collapsed_ncbi_rank_tree, lineage_rank_dict)


    if outpath:
        collapsed_rank_tree.write(outfile=outpath+"_collaped_rank_tree.nwk")
        collapsed_ncbi_rank_tree.write(outfile=outpath+"_collapsed_ncbi_rank_tree.nwk")
        lineage_tree.write(outfile=outpath+"_lineage_tree.nwk")
        ncbi_lineage_tree.write(outfile=outpath + "_ncbi_lineage_tree.nwk")
        lineage_rank_tree.write(outfile=outpath + "_lineage_rank_tree.nwk")
        unique_lineage_rank_tree.write(outfile=outpath + "_unique_lineage_rank_tree.nwk")
        ncbi_lineage_rank_tree.write(outfile=outpath + "_ncbi_lineage_rank_tree.nwk")


# test = ["files/taxon_annotate"]
# test = ["files/test"]
# filepaths = ["files/2U1_complete_sstrimmed_cleaned_phyml", "files/2U1_N_terminal_intact_trimmed_cleaned_phyml", "files/2U1_N_terminal_intact_insertions_removed_trimmed_cleaned_phyml" ]
# for path in test:
#     tax_dict_path = path + "_taxon_dict.obj"
#     compare_tree_with_taxon_tree("family", filepath=path + ".nwk", id=True, filepath_to_load_dict=tax_dict_path)

compare_tree_with_taxon_tree("order", filepath="/Users/gabefoley/Dropbox/PhD/Smaller Projects/GRASP tree/24.nwk", filepath_to_load_dict="../notebooks/files/24.obj", outpath="../notebooks/files/24again")