#! /usr/bin/env python3

"""K-mer propagation post-processing.

Create the main index FASTA file and the new tree, i.e., a minimal tree with k-mer annotations.

Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT
"""

import argparse
import glob
import os
import re
import sys

sys.path.append (os.path.dirname (__file__))
import prophylelib as pro


def create_index_tree():
	pass


def create_fasta(propagation_dir, index_fasta_fn, suffix, verbose=False):
	re_name = re.compile (r'(.*)\.' + suffix.replace (r'.', '\.') + r'$')

	fa_fns = glob.glob ("{}/*".format (propagation_dir))
	fa_fns.sort ()

	with open (index_fasta_fn, "w+") as fa_fo:
		for fn in fa_fns:
			if verbose:
				print ("Processing '{}'".format (fn), file=sys.stderr)

			fn2 = fn.split ("/")[-1]

			m = re_name.match (fn2)

			if not m:
				continue

			node_name = m.group (1)

			with open (fn) as f:
				for x in f:
					if len (x) == 0:
						continue
					if x[0] == ">":
						contig_name = x[1:]
						print (">{}@{}".format (node_name, contig_name), end="", file=fa_fo)
					else:
						print (x, end="", file=fa_fo)


def main():
	parser = argparse.ArgumentParser (
		description='K-mer propagation postprocessing – merging FASTA files and tree minimization.')
	parser.add_argument (
		'dir',
		metavar='<propagation.dir>',
		type=str,
		help='directory with FASTA file'
	)
	parser.add_argument (
		'in_tree_fn',
		type=str,
		metavar='<in.tree.nw>',
		help='input phylogenetic tree',
	)
	parser.add_argument (
		'index_fasta_fn',
		type=str,
		metavar='<index.fa>',
		help='output fast file',
	)
	parser.add_argument (
		'out_tree_fn',
		type=str,
		metavar='<out.tree.nw>',
		help='output phylogenetic tree',
	)
	parser.add_argument (
		'-D',
		dest='nondel',
		action='store_true',
		help='Non-deleting propagation',
	)

	args = parser.parse_args ()

	dir_fn = args.dir
	index_fasta_fn = args.index_fasta_fn

	if args.nondel:
		suffix = "full.fa"
	else:
		suffix = "reduced.fa"

	create_fasta (dir_fn, index_fasta_fn, suffix, verbose=False)


if __name__ == "__main__":
	main ()
