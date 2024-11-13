from ete3 import Tree

tree_file = "tree.nwk"
tree = Tree(tree_file)

print("Original tree:")
print(tree)

camel_virus_id = "MN514967.1"

outgroup = tree.search_nodes(name=camel_virus_id)
if outgroup:
    tree.set_outgroup(outgroup[0])
    print(f"Outgroup '{camel_virus_id}' has been set.")
else:
    print(f"Error: Could not find outgroup with ID '{camel_virus_id}'")

print("Tree with root:")
print(tree)

rooted_tree_file = "camel_rooted_tree.nwk"
tree.write(outfile=rooted_tree_file)

print(f"Rooted tree saved to {rooted_tree_file}")