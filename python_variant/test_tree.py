#! /usr/bin/env python3

from tree_formatter import *
t=read_newick("id_tree_bin.newick",format=10)
#t.render("test.png")

for leaf in t.iter_leaves():
	#print (leaf.name, leaf.named_lineage)
	if hasattr(leaf,"fastapath"):
		print(leaf.fastapath)
	print(leaf.features)
	print(leaf.infasta_seqnum)
	print()
