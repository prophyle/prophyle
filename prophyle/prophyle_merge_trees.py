#! /usr/bin/env python3

"""Merge ProPhyle Newick/NHX trees.

	Author: Karel Brinda <kbrinda@hsph.harvard.edu>

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
		tree_to_add=Tree(x, format=DEFAULT_FORMAT)
		if add_prefixes:
			tree_to_add=add_prefix(tree_to_add, i)
		t.add_child(tree_to_add)

	if verbose:
		print("Writing to '{}'".format(output_tree), file=sys.stderr)

	t.write(
			outfile=output_tree,
			features=[],
			format=DEFAULT_FORMAT,
			format_root_node=True,
		)


def parse_args():
		parser = argparse.ArgumentParser(description='Merge multiple Prophyle trees')

		parser.add_argument('in_tree',
				type=str,
				help='input trees',
				nargs='+',
			)

		parser.add_argument('out_tree',
				type=str,
				help='output trees',
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
