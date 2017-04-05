#! /usr/bin/env python3

import sys
import os
import glob
import argparse

def main():
	parser = argparse.ArgumentParser(description='Merge FASTA files')
	parser.add_argument(
			'dir',
			type=str,
			help='directory with FASTA file'
		)
	parser.add_argument(
			'--nondel',
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
	
	os.chdir(dir_fn)
	fa_fns=glob.glob("*.{}".format(suffix))
	fa_fns.sort()
	for fn in fa_fns:
		if verbose:
			print("Processing '{}'".format(fn), file=sys.stderr)
	
		node=fn.replace("."+suffix,"")
	
		with open(fn) as f:
			for x in f:
				if len(x)==0:
					continue
				if x[0]==">":
					print (">{}@{}".format(node,x[1:]),end="")
				else:
					print(x,end="")

if __name__ == "__main__":
	main()
