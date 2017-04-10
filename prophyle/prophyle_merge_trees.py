#! /usr/bin/env python3

"""Merge ProPhyle Newick/NHX trees.

	Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Todo:
	* Add an option for adding prefixes to node names (e.g., name => tree1:name).
"""

import os
import sys
import argparse

from ete3 import Tree

DEFAULT_FORMAT = 1


def merge_trees(input_trees, output_tree, verbose):
	t = Tree(
			name="merge_root",
		)

	for x in input_trees:
		if verbose:
			print("Loading '{}'".format(x), file=sys.stderr)
		tree_to_add=Tree(x, format=DEFAULT_FORMAT)
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

		parser.add_argument('-v',
				help='verbose',
				dest='verbose',
				action='store_true',
			)

		args = parser.parse_args()
		return args


def main():
	args=parse_args()
	merge_trees(
			input_trees=args.in_tree,
			output_tree=args.out_tree,
			verbose=args.verbose,
		)


if __name__ == "__main__":
	main()
