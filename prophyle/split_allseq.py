#! /usr/bin/env python3

"""Split fasta files containing sequences of different organisms (with different TaxID)
	Used for HMP dataset, since all sequences are distributed in a single fasta

	Author: Simone Pignotti <pignottisimone@gmail.com>

Example:
	$ split_allseq.py <input.fa> <output_dir>

TODO:
	* support infasta_offset and base_len (already included in the trees) for index construction
"""

import os, argparse

def split_fs(output_dir_fn):
	i = 0
	while True:
		i += 1
		fn = output_dir_fn + '/' + str(i) + '.fna'
		yield open(fn, 'w'), i

if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description='Split a fasta file containing multiple sequences in multiple files containing one sequence')
	parser.add_argument('input_file', help='Fasta file')
	parser.add_argument('output_dir', help='Output directory')
	args = parser.parse_args()

	input_f = args.input_file
	output_dir_fn = args.output_dir

	start_seq = '>'
	split_f = split_fs(output_dir_fn)
	outfile = None

	with open(input_f, 'r') as fasta:
		for line in fasta:
			if start_seq not in line:
				outfile.write(line)
			else:
				if outfile:
					outfile.close()
				outfile, i = next(split_f)
				outfile.write(line)

	outfile.close()
