#! /usr/bin/env python3

import sys
import os
import glob

if len(sys.argv)!=2:
	print ("Script has 1 parameter: directory with FASTA files")
	sys.exit(-1)

dir_fn=sys.argv[1]

os.chdir(dir_fn)
fa_fns=glob.glob("*.reduced.fa")
fa_fns.sort()
for fn in fa_fns:
	print("Processing '{}'".format(fn), file=sys.stderr)

	node=fn.replace(".reduced.fa","")

	with open(fn) as f:
		for x in f:
			x=x.strip()
			if len(x)==0:
				continue
			if x[0]==">":
				print (">{}_{}".format(node,x[1:]))
			else:
				print(x)


