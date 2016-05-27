#! /usr/bin/env python3

import sys
import os
import glob

if len(sys.argv)!=2:
    print ("Script has 1 parameter: directory with FASTA files")
    sys.exit(-1)

dir_fn=sys.argv[1]

os.chdir(dir_fn)
for fn in glob.glob("*.reduced.fa"):
    print("Processing '{}'".format(fn), file=sys.stderr)

    parts=fn.split(".")
    contig=parts[0]

    with open(fn) as f:
        for x in f:
            x=x.strip()
            if x=="":
                continue
            if x[0]==">":
                print (">{}_{}".format(contig,x[1:]))
            else:
                print(x)


