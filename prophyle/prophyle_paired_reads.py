#! /usr/bin/env python3

"""Merge paired end fastq files

Author: Simone Pignotti <pignottisimone@gmail.com>

Licence: MIT
"""

import sys
import argparse
import bz2
import gzip
from itertools import islice

desc="""\
	Program: prophyle_paired_reads

	Merge paired end fastq files, eventually compressed with gzip or bz2
"""

magic_dict = {
	b'\x1f\x8b\x08': lambda x: gzip.open(x, 'rt'),
	b'\x42\x5a\x68': lambda x: bz2.open(x, 'rt'),
}

def open_any_file(filename):
	# Initialise to default open function
	max_len = max(len(x) for x in magic_dict)
	with open(filename, 'rb') as f:
		file_start = f.read(max_len)
		f.seek(0)
	for magic, open_func in magic_dict.items():
		if file_start.startswith(magic):
			return open_func(filename)
	return open(filename, 'rt')

def merge_fastas(fastq_1, fastq_2, output_file, log_file):

	i = 1
	while True:
		# get next four lines of each file
		next_read_1 = list(islice(fastq_1, 4))
		next_read_2 = list(islice(fastq_2, 4))
		if not (next_read_1 and next_read_2):
			break
		if len(next_read_1) != 4 or len(next_read_2) != 4:
			print('[prophyle_paired_reads] Error: malformed fastq files ' +
					'(less than 4 lines for read at line ' + str(i) + ')', file=log_file)
			sys.exit(1)
		try:
			id_1 = next_read_1[0].split('/')[0].strip()
			id_2 = next_read_2[0].split('/')[0].strip()
		except KeyError:
			try:
				id_1 = next_read_1[0].split()[0].strip()
				id_2 = next_read_2[0].split()[0].strip()
			except KeyError:
				print('[prophyle_paired_reads] Error: malformed fastq files ' +
						'(cannot split read id at line ' + str(i) + ')', file=log_file)
				sys.exit(1)
		if id_1 != id_2:
			print('[prophyle_paired_reads] Error: malformed fastq files ' +
					'(read ids at line ' + str(i) + ' not corresponding)', file=log_file)
			sys.exit(1)
		out_read = id_1 + '\n' +\
					next_read_1[1].strip() + 'N' + next_read_2[1].strip() + '\n' +\
					next_read_1[2].strip() + '\n' +\
					next_read_1[3].strip() + ':' + next_read_2[3].strip() + '\n'
		output_file.write(out_read)
		i += 4


def main():

	parser = argparse.ArgumentParser(
		description=desc)

	parser.add_argument('fastq_1',
						type = str,
						help = '1st fastq file')
	parser.add_argument('fastq_2',
						type = str,
						help = '2nd fastq file')
	parser.add_argument('-o', '--output-file',
						type=argparse.FileType('w'),
						default = sys.stdout,
						metavar = 'output_file',
						dest = 'output_file',
						help = 'output file [stdout]')
	parser.add_argument('-l', '--log',
						type=argparse.FileType('w'),
						default = sys.stderr,
						metavar = 'log_file',
						dest = 'log_file',
						help = 'log file [stderr]')

	args = parser.parse_args()
	fastq_1 = open_any_file(args.fastq_1)
	fastq_2 = open_any_file(args.fastq_1)
	merge_fastas(fastq_1, fastq_2, args.output_file, args.log_file)
	fastq_1.close()
	fastq_2.close()

if __name__ == '__main__':
	main()
