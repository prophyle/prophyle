#! /usr/bin/env python3

import os
import metag
import datetime

from tree_formatter import *

import logging

logger = logging.getLogger()
handler = logging.StreamHandler()
formatter = logging.Formatter(
        '%(asctime)s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.INFO)

#t=read_newick("id_tree_bin.newick",format=10)
#t.render("test.pdf")

#for leaf in t.iter_leaves():
#	#print (leaf.name, leaf.named_lineage)
#	if hasattr(leaf,"fastapath"):
#		print(leaf.fastapath)
#	print(leaf.features)
#	print(leaf.infasta_seqnum)
#	print()

#def NodeIndex:
#	def __init__(self,node,output_fasta_fn):
#		self.node=node
#		self.output_fasta_fn=output_fasta_fn
#		self.kmers_set=set()
#		self.fasta_fo=None
#
#	def add_kmer(self,kmer):
#		self.kmers_set.add(kmers_set)
#
#	def add_kmers(self,kmers_set):
#		self.kmers_set|=kmers_set
#
#	def __del__(self):
#		if len(self.kmers)>0:
#			metag.set_to_fasta(
#					fasta_fn=self.output_fasta_fn,
#					set_of_kmers=self.kmers_set,
#					assemble=True
#				)


class TreeIndex:

	def __init__(self,tree_newick_fn,format=10,directory="./"):
		self.tree=read_newick(tree_newick_fn,format=format)
		self.directory=directory
		os.makedirs(directory,exist_ok=True)


	#def _rightmost_descendant(self,node):
	#	while not node.is_leaf():
	#		node=node.get_children()[-1]
	#	return node
	#
	#def _rightmost_leaf_of_subtree(self,leaf):
	#	node=leaf
	#	while node.get_ancestors() is not None:
	#		#print(node.name)
	#		ancestors=node.get_ancestors()
	#		if len(ancestors)==0:
	#			break
	#		father=ancestors[0]
	#		fathers_children=father.get_children()
	#		assert(len(fathers_children)<=2)
	#		if len(fathers_children)!=1:
	#			if fathers_children[1]==node:
	#				break
	#		node=father
	#
	#	return _rightmost_descendant(node)
	#
	#def _save_kmers_to_dict(kmers_set,rightmost_leaf):
	#	for kmer in kmers_set:

	@staticmethod
	def _node_debug(node):
		if hasattr(node,"common_name") and node.common_name!="":
			return "{}_{}".format(node.name,node.common_name)
		elif hasattr(node,"sci_name") and node.sci_name!="":
			return "{}_{}".format(node.name,node.sci_name)
		else:
			return "{}".format(node.name)

	def create_fasta(self,node,kmers_set):
		fasta_fn=os.path.join(self.directory,"{}.fa".format(node.name))
		logger.info('Creating FASTA "{}"'.format(fasta_fn))
		logger.debug('... from k-mers "{}"'.format(", ".format(kmers_set)))
		metag.set_to_fasta(
				fasta_fn=fasta_fn,
				set_of_kmers=kmers_set,
				assemble=True,
				contig_prefix="node{}".format(node.name),
			)
		logger.info('... FASTA created "{}"'.format(fasta_fn))

	def get_shared_kmers(self,node,k):
		logger.info('BEGIN get shared k-mers for node "{}"'.format(self._node_debug(node)))
		if node.is_leaf():
			kmers_set=set()
			if hasattr(node,"fastapath"):
				fastas_fn=node.fastapath.split("@")
				#print( datetime.datetime.now().time() )
				for fasta_fn in fastas_fn:
					if os.path.isfile(fasta_fn):
						logger.info(' ...reading "{}"'.format(fasta_fn))
						kmers_set|=metag.set_from_fasta(fasta_fn,k)
			logger.info('END get shared k-mers for node "{}"'.format(self._node_debug(node)))
			logger.debug('... kmers (from fasta files): "{}"'.format(", ".join(kmers_set)))
			return kmers_set
		else:
			children=node.get_children()
			list_of_full_sets=[self.get_shared_kmers(child,k=k) for child in children]
			intersection=set.intersection(*list_of_full_sets)
			list_of_reduced_sets=[full_set-intersection for full_set in list_of_full_sets]
			for (i,reduced_set) in enumerate(list_of_reduced_sets):
				if len(reduced_set)>0:
					self.create_fasta(children[i],reduced_set)
			logger.info('END get shared k-mers for node "{}"'.format(self._node_debug(node)))
			logger.debug('... kmers (from children): "{}"'.format(", ".join(intersection)))
			return intersection

	def build_index(self,k):
		logger.info('Going to build index for k={}'.format(k))
		self.get_shared_kmers(self.tree.get_tree_root(),k=k)


ti=TreeIndex("id_tree_bin.newick",directory="index")
ti.build_index(k=20)
