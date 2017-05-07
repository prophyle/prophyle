#! /usr/bin/env python3

"""Merge paired end fastq files

Author: Simone Pignotti <pignottisimone@gmail.com>

Licence: MIT
"""

import sys
import argparse
from itertools import islice

desc="""\
	Program: prophyle_paired_reads

	Merge paired end fastq files
"""

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
						type = argparse.FileType('r'),
						default = sys.stdin,
						help = 'first read file in the fastq format [stdin]')
	parser.add_argument('fastq_2',
						type = argparse.FileType('r'),
						default = sys.stdin,
						help = 'second read file in the fastq format [stdin]')
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
	merge_fastas(args.fastq_1, args.fastq_2, args.output_file, args.log_file)

if __name__ == '__main__':
	main()
