#! /usr/bin/env python3

import os, sys

names_dmp=sys.argv[1]
ltr_f=sys.argv[2]

records = {}
with open(names_dmp) as dmp:
	for line in dmp:
		values=line.split("|")
		tax_id=int(values[0])
		name=values[1].strip()
		records[tax_id]=name

with open(ltr_f) as ha:
	for line in ha:
		#print(line)
		parts=line.split("\t")
		print(parts[0],end="")
		names=["'{}'".format(records[int(x)]) for x in parts[1:]]
		print("\t".join(names))
		