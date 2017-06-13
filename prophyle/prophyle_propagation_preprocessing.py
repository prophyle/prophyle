#! /usr/bin/env python3

"""K-mer propagation pre-processing.

Extract subtrees, merge ProPhyle Newick/NHX trees and subsample them.

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT

Examples:

		$ prophyle_propagation_preprocessing.py ~/prophyle/bacteria.nw ~/prophyle/viruses.nw bv.nw
		$ prophyle_propagation_preprocessing.py ~/prophyle/bacteria.nw@562 ecoli.nw
"""

import argparse
import ete3
import os
import random
import sys

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro

random.seed(42)


def add_prefix(tree, prefix):
	for node in tree.traverse("postorder"):
		node.name = "{}-{}".format(prefix, node.name)
	return tree


def merge_trees(input_trees_fn, output_tree_fn, verbose, add_prefixes, sampling_rate):
	assert sampling_rate is None or 0.0 <= float(sampling_rate) <= 1.0

	t = ete3.Tree(
		name="merge_root",
	)

	if len(input_trees_fn) == 1:
		if verbose:
			print("Only one tree, don't add any prefix", file=sys.stderr)
		add_prefixes = False

	for i, x in enumerate(input_trees_fn, 1):
		if verbose:
			print("Loading '{}'".format(x), file=sys.stderr)

		tree_fn, _, root_name = x.partition("@")
		tree_to_add = pro.load_nhx_tree(tree_fn)

		# subtree extraction required
		if root_name != '':
			tree_to_add = tree_to_add & root_name

		# prepend prefixes to node names
		if add_prefixes:
			tree_to_add = add_prefix(tree_to_add, i)

		t.add_child(tree_to_add)

	if sampling_rate is not None:
		sampling_rate = float(sampling_rate)

		leaves_1 = []

		for node in t.traverse("postorder"):
			if len(node.children) == 0:
				leaves_1.append(node)

		leaves_2 = random.sample(leaves_1, max(round(sampling_rate * len(leaves_1)), 1))

		leaves_to_remove = set(leaves_1) - set(leaves_2)

		if verbose:
			print("Removing the following leaves: {}".format(", ".join(map(apply(lambda x: x.name, leaves_to_remove)))),
				file=sys.stderr)

		for node in leaves_to_remove:
			while len(node.up.children) == 1:
				node = node.up
			node.detach()

		print("Subsampling the tree with rate {:.4f}, {} leaves out of {} were sampled".format(sampling_rate,
			len(leaves_2), len(leaves_1)), file=sys.stderr)

	if verbose:
		print("Writing to '{}'".format(output_tree), file=sys.stderr)

	pro.save_nhx_tree(t, output_tree_fn)


def parse_args():
	parser = argparse.ArgumentParser(
		description="\n".join(
			[
				'Merge multiple ProPhyle trees. Specific subtrees might be extracted before merging. Examples:',
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

	parser.add_argument('-s',
		help='rate of sampling the tree [no sampling]',
		dest='sampling_rate',
		metavar='FLOAT',
		type=str,
		default=None,
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
	args = parse_args()
	merge_trees(
		input_trees_fn=args.in_tree,
		output_tree_fn=args.out_tree,
		verbose=args.verbose,
		add_prefixes=args.add_prefixes,
		sampling_rate=args.sampling_rate,
	)


if __name__ == "__main__":
	main()
