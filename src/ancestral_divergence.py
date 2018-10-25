def getRemovedSeqs(tree1, tree2):
    removedSeqs = [x.name for x in tree1 if x.name not in tree2]
    removedSeqs2 = [x.name for x in tree2 if x.name not in tree1]

    if len(removedSeqs) > 0 and len(removedSeqs2) > 0:
        raise RuntimeError("There are unique sequences in each tree. One tree should be a subset of the other")

    if len(removedSeqs) > 0:
        return removedSeqs
    elif len(removedSeqs2) > 0:
        return removedSeqs2
    else:
        raise RuntimeError("These trees contained the same sequences")

    return removedSeqs


def getNodeToLeafDict(tree):
    leaves = tree.get_cached_content()
    node2leaf = {}

    for n in tree.traverse():
        if len(leaves[n]) > 1:
            leaf_names = [x.name for x in leaves[n]]
            node2leaf[n] = leaf_names

    return node2leaf


def getLeafToNodeDict(tree):
    leaves = tree.get_cached_content()
    leaf2node = {}

    for n in tree.traverse():
        if len(leaves[n]) > 1:
            leaf_names = [x.name for x in leaves[n]]
            leaf2node[repr(sorted(leaf_names))] = n

    return leaf2node


def mapNodesBetweenTrees(tree1, tree2, removedSeqs):
    nodeDict = {}
    if len(tree1) > len(tree2):
        larger = tree1
        smaller = tree2
    else:
        larger = tree2
        smaller = tree2

    larger_dict = getNodeToLeafDict(larger)
    smaller_dict = getLeafToNodeDict(smaller)

    for node, leaf_names in larger_dict.items():

        removed_names = [x for x in leaf_names if x not in removedSeqs]

        if removed_names in larger_dict.values() and len(removed_names) < len(leaf_names):
            pass
        else:
            leaves = repr(sorted(removed_names))
            print (node.name)
            print (leaf_names)
            nodeDict[node.name] = smaller_dict[leaves].name

    return nodeDict