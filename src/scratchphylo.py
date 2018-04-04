import ete3

print_partitions = False

tree1 = ete3.Tree("/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and bacteria/180312_fifty_percent_identity/Further_exclusions/Comparison_of_trees/2U1_EGcurated_MAFFT_truncated.nwk", format=1)
tree2 = ete3.Tree("/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and bacteria/180312_fifty_percent_identity/Further_exclusions/Comparison_of_trees/RAxML_MAFFT.nwk", format=1)
tree3 = ete3.Tree("/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and bacteria/180312_fifty_percent_identity/Further_exclusions/Comparison_of_trees/2U1_EGcurated_RAxML_MAFFT.nwk", format=1)

tree1_name = "2U1_EGcurated_MAFFT_truncated"
tree2_name = "RAxML_MAFFT"
tree3_name = "2U1_EGcurated_MAFFT"



def checktree(tree1, tree2, tree1_name, tree2_name):

    rf, max_rf, common_leaves, parts_t1, parts_t2, set1, set2 = tree1.robinson_foulds(
        tree2, unrooted_trees=True)
    print("\nRF distance between %s and %s is %s over a total of %s" % (tree1_name, tree2_name, rf, max_rf))
    if print_partitions:
        print("Partitions in %s that were not found in %s: %s" % (tree1_name, tree2_name, parts_t2 - parts_t1))
        print("\nPartitions in %s that were not found in %s: %s" % (tree2_name, tree1_name, parts_t1 - parts_t2))


checktree(tree1,tree2,tree1_name,tree2_name)
checktree(tree1,tree3,tree1_name,tree3_name)
checktree(tree2,tree3,tree2_name,tree3_name)