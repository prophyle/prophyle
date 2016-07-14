#! /usr/bin/env python3

import os
import shutil
import datetime
import sys
import argparse
import operator


from ete3 import Tree

import logging

DEFAULT_FORMAT = 1

class TreeIndex:

	def __init__(self,tree_newick_fn,format=DEFAULT_FORMAT):
		self.tree_newick_fn=tree_newick_fn
		self.tree=Tree(tree_newick_fn,format=1)

		self.name_dict={}

		for node in self.tree.traverse("postorder"):
			self.name_dict[node.name]=node

		print (self.name_dict)


	def assdict_to_weightdic(self,ass_dict):
		all_nodes_hit=set()

		w=ass_dict.copy()

		for node_name in ass_dict:

			node=self.name_dict[node_name]

			while node.up:
				node=node.up

				if node.name in ass_dict:
					w[node_name]+=ass_dict[node.name]

		return w



if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Implementation of assignment algorithm')

	parser.add_argument('-i', '--input',
			type=argparse.FileType('r'),
			required=True,
			dest='input_file',
			help = 'input file'
		)

	parser.add_argument('-n', '--newick-tree',
			type=str,
			metavar='str',
			required=True,
			dest='newick_fn',
			help = 'newick tree',
		)


	args = parser.parse_args()

	newick_fn=args.newick_fn
	inp_fo=args.input_file


	ti=TreeIndex(
			tree_newick_fn=newick_fn,
		)

	print("lets go")

	#ti.process_node(ti.tree.get_tree_root())
	for x in inp_fo:
		stat,rname,ass_node,rlen,ass_kmers=x.split("\t")

		ass_dict={}

		ass_blocks=ass_kmers.split(" ")
		for b in ass_blocks:
			(ids,count)=b.split(":")
			ids=ids.split(",")

			for iid in ids:
				if iid !="0":
					try:
						ass_dict[iid]+=int(count)
					except KeyError:
						ass_dict[iid]=int(count)

		if ass_dict!={}:

			weight_dict=ti.assdict_to_weightdic(ass_dict)
			stat="C"


			#arg max
			#FIX: when 2 arg max exist
			ass_node=max(weight_dict.items(), key=operator.itemgetter(1))[0]

		print(ass_dict, weight_dict)

		print("\t".join([stat,rname,ass_node,rlen,ass_kmers]))
