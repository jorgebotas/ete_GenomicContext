from ete3 import Tree


t = Tree("001_754_949.txt")


for node in t.traverse():
    if node.is_leaf():
        name = node.name
        node.name = name.split(".")[1]


t.write(outfile="001_754_949-tree.txt")


t = Tree("001_754_949-tree.txt")


for node in t.traverse():
    if node.is_leaf():
        name = node.name
        print(name)
