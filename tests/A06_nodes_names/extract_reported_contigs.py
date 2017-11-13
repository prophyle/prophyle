#! /usr/bin/env python3

import fileinput

for line in fileinput.input():
    _, _, _, _, blocks = line.strip().split("\t")
    #print(blocks)
    for x in blocks.split(" "):
        #print(x)
        contigs, _ = x.split(":")
        for y in contigs.split(","):
            print(y)
