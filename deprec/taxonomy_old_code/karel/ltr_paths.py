#! /usr/bin/env python3

import os, sys

taxid_dmp=sys.argv[1]
haha=sys.argv[2]

records = {}
with open(taxid_dmp) as dmp:
	for line in dmp:
		values=line.split("|")
		tax_id=int(values[0])
		parent_tax_id=int(values[1])
		rank=values[2].strip()

		records[tax_id]=(parent_tax_id,rank)

with open(haha) as ha:
	for line in ha:
		#print(line)
		parts=line.split("\t")
		seqname=parts[0].strip()
		taxid=int(parts[2])


		path=[]

		print(seqname,"\t",end="")
		while taxid!=1:
			path.append(taxid)
			taxid=records[taxid][0]
		path.append(taxid)
		print("\t".join(map(str,path[::-1])))
