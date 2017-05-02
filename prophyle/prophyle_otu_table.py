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

	return complete_tree

def create_otu_tables(tree, input_files, target_ranks, field, log):

	# index starts from 0
	taxid_field = field - 1
	# for each node n, tree_ranked_ancestors[n.taxid][r] is its ancestor at rank r
	taxa = [str(n.taxid) for n in tree.traverse('postorder')]
	tree_ranked_ancestors = {t:{} for t in taxa}
	otu_rank_taxacount = {r:Counter() for r in target_ranks}
	otu_tables = {}
	for r in target_ranks:
		for f in input_files:
			try:
				otu_tables[r][f] = Counter()
			except KeyError:
				otu_tables[r] = {}
				otu_tables[r][f] = Counter()

	for node in tree.traverse("postorder"):
		node_taxid = str(node.taxid)
		node_rank = str2rank.get(node.rank,Rank.NO_RANK)
		if node_rank in target_ranks:
			tree_ranked_ancestors[node_taxid][node_rank] = node_taxid
		for anc in node.get_ancestors():
			anc_taxid = str(anc.taxid)
			anc_rank = str2rank.get(anc.rank,Rank.NO_RANK)
			if anc_rank in target_ranks:
				tree_ranked_ancestors[node_taxid][anc_rank] = anc_taxid

	tree_ranked_ancestors['0'] = {}
	tree_ranked_ancestors['*'] = {}
	already_ignored = set()

	for f in input_files:
		with open(f, 'r') as in_f:
			for line in in_f:
				taxid = line.split('\t')[taxid_field].strip()
				try:
					for r,t in tree_ranked_ancestors[taxid].items():
						otu_tables[r][f][t] += 1
						otu_rank_taxacount[r][t] += 1
				except KeyError:
					if taxid not in already_ignored:
						print('[prophyle_otu_table] Error: ignored taxid ' + taxid +
									' (not in the tree)', file=log)
						already_ignored.add(taxid)
					pass

	return otu_tables, otu_rank_taxacount

def write_tables(otu_tables, otu_rank_taxacount, output_prefix, input_files, log):
	rank2str = {v:k for k,v in str2rank.items()}

	out_string = '\t'+'\t'.join(input_files)

	for rank, taxa_count in otu_rank_taxacount.items():
		with open(output_prefix+'_'+rank2str[rank]+'.tsv', 'w') as out_f:
			rank_out_string = out_string
			for taxid, _ in taxa_count.most_common():
				rank_out_string += '\n' + taxid
				for f in input_files:
					rank_out_string += '\t' + str(otu_tables[rank][f][taxid])
			out_f.write(rank_out_string)

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
	parser.add_argument('input_files',
						nargs='+',
						help = 'input files (outputs of prophyle classify)')
	parser.add_argument('-r', '--ranks',
						type = str,
						default = 'species,genus,family,phylum,class,order,kingdom',
						dest = 'target_ranks',
						help = 'comma separated list of ranks to build the OTU table for '
						 		'[species,genus,family,phylum,class,order,kingdom]')
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
	target_ranks = []

	try:
		str_ranks = args.target_ranks.split(',')
	except:
		print('[prophyle_otu_table] Error while parsing ranks: must be a comma' +
					' separated list', file=sys.stderr)
		sys.exit(1)

	target_ranks = [str2rank[r] for r in map(str.strip, str_ranks)]

	complete_tree = build_complete_tree(args.tree, args.log_file)
	otu_tables, otu_rank_taxacount = create_otu_tables(complete_tree, args.input_files,
									target_ranks, args.field, args.log_file)
	write_tables(otu_tables, otu_rank_taxacount, args.output_prefix, args.input_files, args.log_file)

if __name__ == '__main__':
	main()
