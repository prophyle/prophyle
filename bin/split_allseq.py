#! /usr/bin/env python3

import os, argparse

def read_buf(file, buf_size=536870912):
	buf = file.read(buf_size)+file.readline()
	while buf:
		yield buf
		buf = file.read(buf_size)+file.readline()

def split_fs(output_dir_fn):
	i = 0
	while True:
		i += 1
		fn = output_dir_fn + '/' + str(i) + '.fna'
		yield open(fn, 'w'), i

parser = argparse.ArgumentParser(
	description='Split a file containing all the sequences in multiple files')
parser.add_argument(
			'-i','--input-file',
			type=str,
			metavar='str',
			dest='allseq_f',
			required=True,
			help='Fasta file containing all the sequences',
		)
parser.add_argument(
			'-o','--output-dir',
			type=str,
			metavar='str',
			dest='output_dir_fn',
			required=True,
			help='Output directory',
		)
args = parser.parse_args()

allseq_f = args.allseq_f
output_dir_fn = args.output_dir_fn

start_seq = '>'
split_f = split_fs(output_dir_fn)
outfile = None

with open(allseq_f, 'r') as fasta:
	for buf in read_buf(fasta):
		for line in buf:
			if start_seq not in line:
				outfile.write(line)
			else:
				if outfile:
					outfile.close()
				outfile, i = next(split_f)
				outfile.write(line)

outfile.close()
