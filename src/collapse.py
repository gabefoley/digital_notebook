from ete3 import Tree
def collapsed_leaf(node):
    if len(node2labels[node]) == 1:
       return True
    else:
       return False

t = Tree("((((a,a,a)a,a)aa, (b,b)b)ab, (c, (d,d)d)cd);", format=1)
# print (t)
# We create a cache with every node content
node2labels = t.get_cached_content(store_attr="name")
print (node2labels)
# print (t.write(is_leaf_fn=collapsed_leaf))

t2 = Tree( t.write(is_leaf_fn=collapsed_leaf) )
print (t2)