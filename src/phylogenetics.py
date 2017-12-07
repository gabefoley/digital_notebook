import ete3
from ete3 import Tree

latestTree = ete3.Tree("files/2U1_unique_cleaned_C_5.nwk")


# Pull the taxonomy information appended to the end of each ID (thanks SeqScrub!)
taxonomySet = set()
for leaf in latestTree.get_leaf_names():
    taxonomy = leaf.split("_")[-1]
    if taxonomy[0].isupper():
        taxonomySet.add(taxonomy)
print (taxonomySet)
print (latestTree)