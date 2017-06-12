#! /usr/bin/env python3

"""Create the main index FASTA file and move sequences from singletons back down.

Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT
"""

import argparse
import glob
import os
import sys

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro


#def print_singleton_corrections(renaming_dict, file):
#	updates=["'{} => '{}'".format(k, v) if k!=v else "" for (k, v) in renaming_dict.items()]
#	updates=list(filter(lambda x: x!="", updates))
#
#	if len(updates)!=0:
#		print("Singleton corrections:", ", ".join(updates), file=file)


#def build_renaming_dict(tree_fn):
#	ren_dict={}
#	tree=pro.load_nhx_tree(tree_fn, validate=True)
#
#	for node in tree.traverse("postorder"):
#		orig_node=node
#		new_node=node
#		while len(new_node.children)==1:
#			new_node = new_node.children[0]
#		ren_dict[orig_node.name]=new_node.name
#	print_singleton_corrections(ren_dict, file=sys.stderr)
#	return ren_dict


def main():
	parser = argparse.ArgumentParser(description='Merge FASTA files')
	parser.add_argument(
			'dir',
			type=str,
			metavar='<propagation.dir>',
			help='directory with FASTA file'
		)
	parser.add_argument(
			'newick_fn',
			type=str,
			metavar='<tree.nw>',
			help='phylogenetic tree for singleton correction (in Newick/NHX) [None]',
			default = None,
		)
	parser.add_argument(
			'-D',
			dest='nondel',
			action='store_true',
			help='Non-deleting propagation',
		)

	verbose=False
	
	args = parser.parse_args()
	
	dir_fn=args.dir
	if args.nondel:
		suffix="full.fa"
	else:
		suffix="reduced.fa"
	
	#if args.newick_fn is None:
	#	node_rename_func = lambda x: x
	#else:
	#	ren_dict=build_renaming_dict(args.newick_fn)
	#	node_rename_func = lambda x: ren_dict[x]

	os.chdir(dir_fn)
	fa_fns=glob.glob("*.{}".format(suffix))
	fa_fns.sort()
	for fn in fa_fns:
		if verbose:
			print("Processing '{}'".format(fn), file=sys.stderr)
	
		node_name=fn.replace("."+suffix,"")
		#node_name=node_rename_func(node_name)
	
		with open(fn) as f:
			for x in f:
				if len(x)==0:
					continue
				if x[0]==">":
					contig_name=x[1:]
					print (">{}@{}".format(node_name, contig_name), end="")
				else:
					print (x, end="")


if __name__ == "__main__":
	main()
