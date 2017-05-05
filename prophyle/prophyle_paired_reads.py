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

    while True:
        try:
            next_n_lines_1 = list(islice(fastq_1, 4))
            next_n_lines_2 = list(islice(fastq_2, 4))
            if not (next_n_lines_1 and next_n_lines_2):
                break
            else:
                assert len(next_n_lines_1) == 4 and len(next_n_lines_2) == 4
				try:
					id_1 = next_n_lines_1[0].split('/')[0]
					id_2 = next_n_lines_2[0].split('/')[0]
				except:
					try:
						id_1 = next_n_lines_1[0].split()[0]
						id_2 = next_n_lines_2[0].split()[0]
					except:
						print('[prophyle_paired_reads] Error: malformed fastq files or paired end reads non matching')
						sys.exit(1)
				assert id_1 == id_2
				out_read = id_1 + '\n' + next_n_lines_1[1] + 'N' + next_n_lines_2[1]
					+ '\n' + next_n_lines_1[2] + '\n' + next_n_lines_1[3] + ':' + next_n_lines_2[3]
				output_file.write(out_read)
        except AssertionError:
            print('[prophyle_paired_reads] Error: malformed fastq files or paired end reads non matching')
			sys.exit(1)
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
