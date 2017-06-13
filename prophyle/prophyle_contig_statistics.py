#! /usr/bin/env python3

import argparse
import os
import statistics
import sys

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro


def main():
	parser = argparse.ArgumentParser(
		description='Compute contig statistics for a ProPhyle index. Fasta index should be already computed (using samtools faidx).')
	parser.add_argument(
		'-k',
		type=int,
		metavar='int',
		dest='k',
		default=None,
		help='k-mer length [detected automatically]',
	)
	parser.add_argument(
		'index_dir',
		metavar='<index.dir>',
		type=str,
		help='index directory (will be created)',
	)

	args = parser.parse_args()

	lengths = []

	fai_fn = os.path.join(args.index_dir, "index.fa.fai")

	if args.k is None:
		k = pro.detect_k_from_index(args.index_dir)
	else:
		k = args.k

	with open(fai_fn) as file:
		for x in file:
			_, length, _, _, _ = x.split()
			lengths.append(int(length))

		lengths.sort()

		contig_nb = len(lengths)
		len_total = sum(lengths)
		len_mean = statistics.mean(lengths)
		len_stdev = statistics.stdev(lengths)
		len_median = statistics.median(lengths)

		kmer_occ = len_total - contig_nb * (k - 1)

		print("Number of contigs:\t{}".format(contig_nb))
		print("Total length:\t{}".format(len_total))
		print("Average length:\t{:.2f}".format(len_mean))
		print("Standard devation:\t{:.2f}".format(len_stdev))
		print("Median length:\t{}".format(len_median))
		print("Number of k-mer occurencies:\t{}".format(kmer_occ))


if __name__ == "__main__":
	main()
