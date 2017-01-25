#! /usr/bin/env python3

import os
import re
import sys
import argparse

from subprocess import Popen, PIPE

from ete3 import Tree

#FEATURES=["lineage", "named_lineage", "seqname", "dist", "name",
#	"support", "taxid", "rank", "base_len", "fastapath",
#	"sci_name", "infasta_offset", "gi",
#	'kmers_full','kmers_reduced']


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
					assert counts[cat][nname]==count, "Different k-mer sizes reported for the same node '{}', probably a bug of prophyle-assembler".format(nname)
				except KeyError:
					counts[cat][nname]=count
	return counts

def enrich_tree(
		inp_tree_fn,
		out_tree_fn,
		count_tb,
	):


	tree=Tree(inp_tree_fn,format=1)

	# compute k-mer counts
	#print(sorted(count_tb["full"].keys()))
	#print(sorted(count_tb["reduced"].keys()))
	for node in tree.traverse("preorder"):
		node.del_feature('kmers_full')
		node.del_feature('kmers_reduced')

		nname=node.name

		# todo: nodes with name="" should not exist

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

	# regularly update
	tree.write(
			format=1,
			features=[],
			outfile=out_tree_fn,
			format_root_node=True,
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
	#parser.add_argument(
	#		'-k',
	#		type=int,
	#		metavar='int',
	#		dest='k',
	#		required=True,
	#		help='K-mer length k.',
	#	)
	parser.add_argument(
			'-c','--counts',
			type=str,
			metavar='str',
			dest='counts_fn',
			required=True,
			help='TSV file with counts produced during index propagation.',
		)

	args = parser.parse_args()

	#k=args.k
	#assert k>0
	inp_tree_fn=args.inp_tree_fn
	out_tree_fn=args.out_tree_fn
	counts_fn=args.counts_fn
	
	count_tb=load_nb_kmers(counts_fn)

	enrich_tree(
		#k=k,
		inp_tree_fn=inp_tree_fn,
		out_tree_fn=out_tree_fn,
		count_tb=count_tb,
	)
	
