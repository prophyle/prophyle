#! /usr/bin/env python3

import os, argparse

def split_fs(output_dir_fn):
	i = 0
	while True:
		i += 1
		fn = output_dir_fn + '/' + str(i) + '.fna'
		yield open(fn, 'w'), i

parser = argparse.ArgumentParser(
	description='Split a file containing all the sequences in multiple files')
parser.add_argument('input_file', help='Fasta file containing all the sequences')
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
