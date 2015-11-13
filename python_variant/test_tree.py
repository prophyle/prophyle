#! /usr/bin/env python3

import metag
from tree_formatter import *
t=read_newick("id_tree_bin.newick",format=10)
#t.render("test.pdf")

#for leaf in t.iter_leaves():
#	#print (leaf.name, leaf.named_lineage)
#	if hasattr(leaf,"fastapath"):
#		print(leaf.fastapath)
#	print(leaf.features)
#	print(leaf.infasta_seqnum)
#	print()

def find_rightmost_leaf(node):
	while not node.is_leaf():
		node=node.get_children()[-1]
	return node

def find_last_possible_leaf(leaf):
	node=leaf
	while node.get_ancestors() is not None:
		#print(node.name)
		ancestors=node.get_ancestors()
		if len(ancestors)==0:
			break
		father=ancestors[0]
		fathers_children=father.get_children()
		assert(len(fathers_children)<=2)
		if len(fathers_children)!=1:
			if fathers_children[1]==node:
				break
		node=father
	return find_rightmost_leaf(node)

for leaf in t.iter_leaves():
	last_possible=find_last_possible_leaf(leaf)
	print(leaf.name,last_possible.name)
	if hasattr(leaf,"fastapath"):
		print(leaf.fastapath)
		fastas_fn=leaf.fastapath.split("@")
		kmers_set=set()
		for fasta_fn in fastas_fn:
			kmers_set|=metag.set_from_fasta(fasta_fn,10)
		print(kmers_set)

