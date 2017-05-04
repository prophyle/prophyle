#! /usr/bin/env python3

"""Incorporate information about k-mer counts at individual nodes into a ProPhyle NHX tree.

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT
"""


import os
import re
import sys
import argparse

from ete3 import Tree


def load_nb_kmers(tsv_fn):
	re_fa=re.compile(r'(.*)\.([a-z]*)\.fa')

	counts = {"reduced":{}, "full": {}}
	with open(tsv_fn) as f:
		for x in f:
			if len(x)==0:
				continue
			if x[0]=="#":
				continue
			else:
				x = x.strip()
				(fa,count)=x.split("\t")
				count=int(count)
				fa_short = os.path.basename(fa)
				m=re_fa.match(fa_short)

				cat=m.group(2)
				nname=m.group(1)

				try:
					assert counts[cat][nname]==count, "Different k-mer sizes reported for the same node '{}', probably a bug of prophyle_assembler".format(nname)
				except KeyError:
					counts[cat][nname]=count

	return counts


def enrich_tree(
		inp_tree_fn,
		out_tree_fn,
		count_tb,
	):

	tree=Tree(inp_tree_fn,format=1)

	for node in tree.traverse("preorder"):
		node.del_feature('kmers_full')
		node.del_feature('kmers_reduced')

		nname=node.name

		assert nname != "", "There is a node without any name ('')"

		try:
			node.add_features(kmers_full=count_tb["full"][nname])
		except KeyError:
			print("Warning: full-{} is missing".format(nname),file=sys.stderr)
			node.add_features(kmers_full=0)

		try:
			node.add_features(kmers_reduced=count_tb["reduced"][nname])
		except KeyError:
			print("Warning: reduced-{} is missing".format(nname),file=sys.stderr)
			node.add_features(kmers_reduced=0)


		assert 0<=node.kmers_reduced
		assert node.kmers_reduced<=node.kmers_full


	# make saving newick reproducible
	features=set()
	for n in tree.traverse():
		features|=n.features

	# otherwise some names stored twice – also as a special attribute
	features.remove("name")

	# regularly update
	tree.write(
			format=1,
			features=sorted(features),
			outfile=out_tree_fn,
			format_root_node=True,
		)


def main():
	parser = argparse.ArgumentParser(description='Build index.')
	parser.add_argument('tree',
			metavar='<in-tree.nw>',
			type=str,
			help='input phylogenetic tree (in Newick/NHX)',
		)
	parser.add_argument('counts_fn',
			metavar='<counts.tsv>',
			type=str,
			help='TSV file with counts produced during index propagation.',
		)
	parser.add_argument('out_tree_fn',
			metavar='<out-tree.nw>',
			type=str,
			help='output phylogenetic tree (in Newick/NHX)',
		)

	args = parser.parse_args()

	inp_tree_fn=args.inp_tree_fn
	out_tree_fn=args.out_tree_fn
	counts_fn=args.counts_fn
	
	count_tb=load_nb_kmers(counts_fn)

	enrich_tree(
		inp_tree_fn=inp_tree_fn,
		out_tree_fn=out_tree_fn,
		count_tb=count_tb,
	)


if __name__ == "__main__":
	main()