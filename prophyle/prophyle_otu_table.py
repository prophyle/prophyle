#! /usr/bin/env python3

"""Create an OTU table for ProPhyle's classification output

Author: Simone Pignotti <pignottisimone@gmail.com>

Licence: MIT
"""

import sys
import os
import argparse
import enum
from collections import Counter

#EXTERNAL
from ete3 import Tree, NCBITaxa

desc="""\
	Program: prophyle_otu_table

	Create an OTU table for ProPhyle's classification output
"""

class OrderedEnum(enum.Enum):
	def __ge__(self, other):
		if self.__class__ is other.__class__:
			return self.value >= other.value
		return NotImplemented
	def __gt__(self, other):
		if self.__class__ is other.__class__:
			return self.value > other.value
		return NotImplemented
	def __le__(self, other):
		if self.__class__ is other.__class__:
			return self.value <= other.value
		return NotImplemented
	def __lt__(self, other):
		if self.__class__ is other.__class__:
			return self.value < other.value
		return NotImplemented

class Rank(OrderedEnum):
	NO_RANK = 7
	SPECIES = 6
	GENUS = 5
	FAMILY = 4
	ORDER = 3
	CLASS = 2
	PHYLUM = 1
	KINGDOM = 0

str2rank = {
		'species':Rank.SPECIES,
		'genus':Rank.GENUS,
		'family':Rank.FAMILY,
		'order':Rank.ORDER,
		'class':Rank.CLASS,
		'phylum':Rank.PHYLUM,
		'kingdom':Rank.KINGDOM,
		'no rank':Rank.NO_RANK
		}

def build_complete_tree(tree, log):
	"""Build the taxonomic tree including internal nodes (for rank dependent evaluation)

	Args:
		tree (fileobject): file name of the reference tree, without internal nodes
		log (fileobject): log file
	"""
	ncbi = NCBITaxa()
	original_tree = Tree(tree, format=1)
	taxa = [l.taxid for l in original_tree]
	built = False
	while not built:
		try:
			complete_tree = ncbi.get_topology(taxa, intermediate_nodes=True)
			built = True
		except KeyError as e:
			taxid_not_found = int(e.args[0])
			taxa.remove(taxid_not_found)
			if log:
				print('[prophyle_otu_table] ERROR: TaxID ' + str(taxid_not_found) +
						' not found in ETE DB (try updating it)', file=log)
			pass
	# for node in complete_tree.traverse('postorder'):
	# 	node.name = node.taxid
	return complete_tree

def create_otu_tables(tree, input_file, target_ranks, field, log):

	# index starts from 0
	taxid_field = field - 1
	# for each node n, tree_ranked_ancestors[n.taxid][r] is its ancestor at rank r
	taxa = [str(n.taxid) for n in tree.traverse('postorder')]
	tree_ranked_ancestors = {t:{} for t in taxa}
	tree_ranks = {}

	otu_tables = {r:Counter() for r in target_ranks}

	for node in tree.traverse("postorder"):
		node_taxid = str(node.taxid)
		node_rank = str2rank.get(node.rank,Rank.NO_RANK)
		tree_ranks[node_taxid] = node_rank
		if node_rank in target_ranks:
			tree_ranked_ancestors[node_taxid][node_rank] = node_taxid
		for anc in node.get_ancestors():
			anc_taxid = str(anc.taxid)
			anc_rank = str2rank.get(anc.rank,Rank.NO_RANK)
			if anc_rank in target_ranks:
				tree_ranked_ancestors[node_taxid][anc_rank] = anc_taxid

	tree_ranks['0'] = Rank.NO_RANK
	tree_ranked_ancestors['0'] = {}
	tree_ranks['*'] = Rank.NO_RANK
	tree_ranked_ancestors['*'] = {}

	for line in input_file:
		taxid = line.split('\t')[taxid_field].strip()
		rank = tree_ranks[taxid]
		try:
			for r,t in tree_ranked_ancestors[taxid].items():
				otu_tables[r][t] += 1
		except KeyError:
			print('[prophyle_otu_table] Taxid ' + taxid +
						' is not in the tree', file=log)
			pass

	return otu_tables

def write_tables(otu_tables, output_prefix, log):
	rank2str = {v:k for k,v in str2rank.items()}
	for rank, counter in otu_tables.items():
		with open(output_prefix+'_'+rank2str[rank]+'.tsv', 'w') as out_f:
			for taxid, count in counter.most_common():
				print(taxid+'\t'+str(count),file=out_f)

def main():

	parser = argparse.ArgumentParser(
		description=desc)

	parser.add_argument('tree',
						type = str,
						metavar = '<tree>',
						help = 'taxonomic tree used for classification')
	parser.add_argument('output_prefix',
						type = str,
						metavar = '<output_prefix>',
						help = 'prefix for output files (one per rank, each with suffix \"_rank.tsv\")')
	parser.add_argument('input_file',
						nargs='?',
						type = argparse.FileType('r'),
						default = sys.stdin,
						help = 'input file (output of prophyle classify) [stdin]')
	parser.add_argument('-r', '--ranks',
						type = str,
						default = 'species,genus,family,phylum,class,order,kingdom',
						dest = 'target_ranks',
						help = 'ranks to build the OTU table for [species,genus,family,phylum,class,order,kingdom]')
	parser.add_argument('-f','--field',
						type = int,
						default = 3,
						dest = 'field',
						help = 'position of the taxid in the input lines [3]')
	parser.add_argument('-l', '--log',
						type=argparse.FileType('w'),
						default = sys.stderr,
						metavar = 'log_file',
						dest = 'log_file',
						help = 'log file [stderr]')

	args = parser.parse_args()
	unknown_ranks = set()
	target_ranks = []

	try:
		str_ranks = args.target_ranks.split(',')
	except:
		print('[prophyle_otu_table] Error while parsing ranks: must be a comma' +
					' separated list', file=sys.stderr)
		sys.exit(1)

	target_ranks = [str2rank[r] for r in map(str.strip, str_ranks)]

	complete_tree = build_complete_tree(args.tree, args.log_file)
	otu_tables = create_otu_tables(complete_tree, args.input_file,
									target_ranks, args.field, args.log_file)
	write_tables(otu_tables, args.output_prefix, args.log_file)

if __name__ == '__main__':
	main()
