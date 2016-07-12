#! /usr/bin/env python3

import ete3
import argparse
from tree_formatter import *

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Input Newick tree.')
	parser.add_argument(
			'-i','--input-newick-tree',
			type=str,
			metavar='str',
			dest='newick_i_fn',
			required=True,
			help='Input Newick tree.',
		)
	parser.add_argument(
			'-o','--output-newick-tree',
			type=str,
			metavar='str',
			dest='newick_o_fn',
			required=True,
			help='Output Newick tree.',
		)

	args = parser.parse_args()

	t=read_newick(args.newick_i_fn,format=10)
	t.write(outfile=args.newick_o_fn,features=["fastapath","common_name","sci_name","named_lineage","name","infasta_seqnum","est_nkmers","taxid","rank","base_len","seqname"])
