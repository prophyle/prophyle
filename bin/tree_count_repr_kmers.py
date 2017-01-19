#! /usr/bin/env python3

"""
	Parameters:
		- NONPROP: no k-mer propagation (sequences for leaves only)
		- REASM: re-assemble sequences in leaves
		- NONDEL: non-deletative propagation, implies REASM
		- MASKREP: mask repeats in leaves
"""

import os
import shutil
import datetime
import sys
import argparse

from subprocess import Popen, PIPE

from ete3 import Tree


SCRIPT_DIR=os.path.dirname(os.path.realpath(__file__))
ASM=os.path.join(SCRIPT_DIR,"prophyle-assembler")


def get_nb_kmers(fa_fn):
	fh = open("NUL","w")
	process = Popen([ASM, "-i", fa_fn], stdout=PIPE, stderr=fh)
	(output, err) = process.communicate()
	exit_code = process.wait()
	assert exit_code==0
	out=output.decode('ascii').strip()
	parts=out.split("\t")
	#print(parts)
	assert(len(parts)==2)
	return int(parts[1])

def enrich_tree(
		k,
		inp_tree_fn,
		out_tree_fn,
		index_dir,
	):


	tree=Tree(inp_tree_fn,format=1)

	# remove k-mer counts
	for node in tree.traverse("postorder"):
		node.del_feature('no_repr_kmers_full')
		node.del_feature('no_repr_kmers_reduced')


	# compute k-mer counts
	for node in tree.traverse("postorder"):
		nname=node.name
		print (nname)

		if nname != "":
			fn_red=os.path.join(index_dir, "{}.reduced.fa".format(nname))
			fn_full=os.path.join(index_dir, "{}.full.fa".format(nname))

			nb_red=get_nb_kmers(fn_red)
			nb_full=get_nb_kmers(fn_full)

			node.add_features(no_repr_kmers_full=nb_full)
			node.add_features(no_repr_kmers_reduced=nb_red)

	tree.write(
			format=1,
			features=['fastapath','no_repr_kmers_full','no_repr_kmers_reduced'],
			outfile=out_tree_fn,
		)

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Build index.')
	parser.add_argument(
			'-i','--inp-tree',
			type=str,
			metavar='str',
			dest='inp_tree_fn',
			required=True,
			help='Input taxonomic tree (Newick).',
		)
	parser.add_argument(
			'-o','--out-tree',
			type=str,
			metavar='str',
			dest='out_tree_fn',
			required=True,
			help='Output taxonomic tree (Newick).',
		)
	parser.add_argument(
			'-k',
			type=int,
			metavar='int',
			dest='k',
			required=True,
			help='K-mer length k.',
		)
	parser.add_argument(
			'-d','--index-dir',
			type=str,
			metavar='str',
			dest='index_dir',
			required=True,
			help='Index directory (with *.{full,reduced}.fa).',
		)

	args = parser.parse_args()

	k=args.k
	assert k>0
	inp_tree_fn=args.inp_tree_fn
	out_tree_fn=args.out_tree_fn
	index_dir=args.index_dir

	

	enrich_tree(
		k=k,
		inp_tree_fn=inp_tree_fn,
		out_tree_fn=out_tree_fn,
		index_dir=index_dir,
	)
