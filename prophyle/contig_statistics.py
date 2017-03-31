#! /usr/bin/env python3

import argparse
import statistics

parser = argparse.ArgumentParser(description='Input Newick tree.')
parser.add_argument(
		'-k',
		type=int,
		metavar='int',
		dest='k',
		required=True,
		help='k-mer length',
	)
parser.add_argument(
		'-f','--fai',
		type=str,
		metavar='str',
		dest='fai_fn',
		required=True,
		help='Fasta index (.fai).',
	)

args = parser.parse_args()

lengths=[]

with open(args.fai_fn) as file:
	for x in file:
		c = x.split()
		lengths.append(int(c[1]))

	lengths.sort()

	contig_nb=len(lengths)
	len_total=sum(lengths)
	len_mean=statistics.mean(lengths)
	len_stdev=statistics.stdev(lengths)
	len_median=statistics.median(lengths)

	kmer_occ=len_total - contig_nb * (args.k - 1)

	print("Number of contigs: {}".format(contig_nb))
	print("Total length: {}".format(len_total))
	print("Average length: {}".format(len_mean))
	print("   ..st. dev: {}".format(len_stdev))
	print("Median length: {}".format(len_median))
	print("Number of k-mer occurencies: {}".format(kmer_occ))
