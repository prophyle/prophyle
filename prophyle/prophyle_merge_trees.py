#! /usr/bin/env python3

"""Merge ProPhyle Newick/NHX trees.

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT

Examples:

		$ prophyle_merge_trees.py ~/prophyle/bacteria.nw ~/prophyle/viruses.nw bv.nw
		$ prophyle_merge_trees.py ~/prophyle/bacteria.nw@562 ecoli.nw
"""

import os
import sys
import argparse

from ete3 import Tree

DEFAULT_FORMAT = 1


def add_prefix(tree, prefix):
	for node in tree.traverse("postorder"):
		node.name="{}-{}".format(prefix, node.name)
	return tree


def merge_trees(input_trees, output_tree, verbose, add_prefixes):
	t = Tree(
			name="merge_root",
		)

	if len(input_trees)==1:
		add_prefixes=False

	for i,x in enumerate(input_trees,1):
		if verbose:
			print("Loading '{}'".format(x), file=sys.stderr)
		tree_fn,_,root_name=x.partition("@")
		tree_to_add=Tree(tree_fn, format=DEFAULT_FORMAT)
		if root_name!='':
			tree_to_add=tree_to_add&root_name
		if add_prefixes:
			tree_to_add=add_prefix(tree_to_add, i)
		t.add_child(tree_to_add)

	if verbose:
		print("Writing to '{}'".format(output_tree), file=sys.stderr)

	# make saving newick reproducible
	features=set()
	for n in t.traverse():
		features|=n.features

	# otherwise some names stored twice – also as a special attribute
	features.remove("name")

	t.write(
			outfile=output_tree,
			features=sorted(features),
			format=DEFAULT_FORMAT,
			format_root_node=True,
		)


def parse_args():
		parser = argparse.ArgumentParser(
				description="\n".join(
					[
						'Merge multiple Prophyle trees. Specific subtrees might be extracted before merging. Examples:',
						'\t$ prophyle_merge_trees.py ~/prophyle/bacteria.nw ~/prophyle/viruses.nw bv.nw',
						'\t$ prophyle_merge_trees.py ~/prophyle/bacteria.nw@562 ecoli.nw'
					]),
				formatter_class=argparse.RawTextHelpFormatter,
			)

		parser.add_argument('in_tree',
				metavar='<in_tree.nw{@node_name}>',
				type=str,
				help='input tree',
				nargs='+',
			)

		parser.add_argument('out_tree',
				metavar='<out_tree.nw>',
				type=str,
				help='output tree',
			)

		parser.add_argument('-V',
				help='verbose',
				dest='verbose',
				action='store_true',
			)

		parser.add_argument('-P',
				help='do not add prefixes to node names',
				dest='add_prefixes',
				action='store_false',
			)

		args = parser.parse_args()
		return args


def main():
	args=parse_args()
	merge_trees(
			input_trees=args.in_tree,
			output_tree=args.out_tree,
			verbose=args.verbose,
			add_prefixes=args.add_prefixes,
		)


if __name__ == "__main__":
	main()
