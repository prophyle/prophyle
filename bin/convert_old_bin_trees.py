#! /usr/bin/env python3

import sys
from ete3 import Tree

t = Tree(sys.argv[1], format=1)

root = t.children[0]
count = 0
for node in root.traverse("postorder"):
	if hasattr(node, "taxid"):
		count += 1
digits = len(str(count+1))
root.name = "n"+("0"*(digits-len(str(count+1))))+str(count+1)
root.add_features(taxid = "0")

new_id = 1
for node in root.traverse("postorder"):
	if hasattr(node, "taxid"):
		node.name = "n"+("0"*(digits-len(str(new_id))))+str(new_id)
		new_id += 1

for node in root.traverse("postorder"):
	if not hasattr(node,"taxid"):
		fake_depth = 1
		temp = node.up
		while not hasattr(temp,"taxid"):
			temp = temp.up
			fake_depth += 1
		while temp.search_nodes(name=temp.name+"."+str(fake_depth)):
			fake_depth += 1
		node.name = temp.name+"."+str(fake_depth)
		node.add_features(taxid = temp.taxid)
		if hasattr(temp,"lineage"):
			node.add_features(lineage = temp.lineage, named_lineage = temp.named_lineage,
					rank = temp.rank, sci_name = temp.sci_name)

# for node in root.traverse("postorder"):
# 	assert node.name.startswith('n')
# assert root.name.startswith('n')

# for node in root.traverse("postorder"):
# 	if len(root.search_nodes(name=node.name)) > 1:
# 		print("DUPLICATE: " + node.name)

root.write(features = ["lineage", "named_lineage", "seqname", "dist", "name",
					"support", "taxid", "rank", "base_len", "fastapath",
					"sci_name", "infasta_offset", "gi"],
			format = 1,
			format_root_node = True,
			outfile = sys.argv[2])
