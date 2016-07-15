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

		#print (self.name_dict)

	def dict_from_list(self,kmers_assigned_l):
		d={}
		for (noden_l, count) in kmers_assigned_l:
			for noden in noden_l:
				try:
					d[noden]+=count
				except KeyError:
					d[noden]=count
		return d

	def dict_from_list_lca(self,kmers_assigned_l):
		d={}
		for (noden_l, count) in kmers_assigned_l:
			noden=self.lca(noden_l)
			try:
				d[lca.name]+=count
			except KeyError:
				d[lca.name]=count
		return d

	def lca(self,noden_l):
		nodes_l=list(map(lambda x:self.name_dict[x],noden_l))
		lca=nodes_l[0]
		for node in nodes_l:
			lca=lca.get_common_ancestor(node)
		return lca.name

	def assign(self,kmers_assigned_l,simulate_lca=False):
		all_nodes_hit=set()

		if simulate_lca:
			d=self.dict_from_list_lca(kmers_assigned_l)
		else:
			d=self.dict_from_list(kmers_assigned_l)
		w=d.copy()

		for noden in d:
			node=self.name_dict[noden]
			while node.up:
				node=node.up
				if node.name in d:
					w[noden]+=d[node.name]
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

	parser.add_argument('-l', '--sim-lca',
			action='store_true',
			dest='lca',
			help = 'simulate LCA',
		)


	args = parser.parse_args()

	newick_fn=args.newick_fn
	inp_fo=args.input_file
	lca=args.lca


	ti=TreeIndex(
			tree_newick_fn=newick_fn,
		)

	#print("lets go")

	#ti.process_node(ti.tree.get_tree_root())
	for x in inp_fo:
		x=x.strip()
		stat,rname,ass_node,rlen,a_kmers=x.split("\t")

		l=[]

		blocks=a_kmers.split(" ")
		for b in blocks:
			(ids,count)=b.split(":")

			if ids=="A" or ids=="0":
				continue

			l.append((ids.split(","),int(count)))

		if l!=[]:

			a=ti.assign(l,simulate_lca=lca)
			stat="C"

			max_hit=-1
			noden_m_l=[]
			for noden in a:
				if a[noden]==max_hit:
					noden_m_l.append(noden)
				elif a[noden]>max_hit
					noden_m_l=[noden]
					max_hit=a[noden]

		print("\t".join([stat,rname,ass_node,rlen,a_kmers]))
