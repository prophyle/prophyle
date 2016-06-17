#! /usr/bin/env python3

import sys,os

fai_f = sys.argv[1]
taxid_map_f = sys.argv[2]


gis1=[]

with open(fai_f) as f1:
    for line1 in f1:
        (seqname,_,_)=line1.partition("\t")
        values=seqname.split("|")
        gi1=values[2]
        try:
            gis1.append([seqname,int(gi1)])
        except:
            pass

gis1=sorted(gis1,key=lambda x:x[1])
#print(gis1)


with open(taxid_map_f) as taxid_map:
    for line2 in taxid_map:
        (gi2, _, taxid2) = line2.partition("\t")
        gi2=int(gi2)
        while gis1[0][1] < gi2:
            del gis1[0]
            if len(gis1)==0:
                sys.exit(0)
        if gis1[0][1]==gi2:
            print("\t".join(map(str,[gis1[0][0],gis1[0][1],taxid2.strip()])))
