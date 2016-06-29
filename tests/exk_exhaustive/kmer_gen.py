#! /usr/bin/env python3

import sys, os, argparse
from itertools import product, repeat

parser = argparse.ArgumentParser(description =
			'Generate all possible k-mers with given lenght')

parser.add_argument('k', type=int, help = 'k-mers lenght')
parser.add_argument('-o', default = 'gen_kmers.fq',
					help = 'output file (default: gen_kmers.fq)')

args = parser.parse_args()

with open(args.o, 'w') as output:
	i = 1
	for kmer in product('ACGT', repeat=args.k):
		output.write('@seq'+str(i)+'\n'+
					''.join(kmer)+'\n'+
					'+\n'+''.join(repeat('~',args.k))+'\n')
		i += 1
