#! /usr/bin/env python3

"""Merge paired-end FASTA or FASTQ files

Author: Simone Pignotti <pignottisimone@gmail.com>

Licence: MIT
"""

import sys
import os
import argparse

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro


def read_id(read_1, read_2):
	# two possibilities:
	# <instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>   <pair idx>:<is filtered>:<control number>:<index sequence>
	# <instrument>:<flowcell lane>:<tile>:<x-pos>:<y-pos>:<multiplex idx>/<pair idx>
	try:
		id1 = '/'.join(read_1.split('/')[:-1])
		id2 = '/'.join(read_2.split('/')[:-1])
		if id1 != id2:
			raise ValueError('different id')
		else:
			return id1
	except (KeyError, ValueError):
		try:
			id1 = ' '.join(read_1.split(' ')[:-1])
			id2 = ' '.join(read_2.split(' ')[:-1])
			if id1 != id2:
				raise ValueError('different id')
			else:
				return id1
		except (KeyError, ValueError):
			# cannot find id, just use the one of the first read
			return read_1


def merge_reads(f_reads_1, f_reads_2):

	reads_1 = pro.open_gzip(f_reads_1)
	reads_2 = pro.open_gzip(f_reads_2)

	first_read_1 = reads_1.readline()
	first_read_2 = reads_2.readline()

	first_char_1 = first_read_1[0]
	first_char_2 = first_read_2[0]
	assert first_char_1 == first_char_2, 'paired-end files of different format'
	if first_char_1 == '@':
		fastq, fasta = True, False
	elif first_char_1 == '>':
		fastq, fasta = False, True
	else:
		print('[prophyle_paired_reads] Error: unknown read format (does not start with > or @)',
				file=sys.stderr)
		sys.exit(1)

	print(read_id(first_read_1,first_read_2))

	if fastq:
		# file line
		i = 0
		for next_line_1, next_line_2 in zip(reads_1, reads_2):
			i += 1
			l = i%4
			if l == 0:
				assert next_line_1.startswith('@') and next_line_2.startswith('@'), \
					"malformed fastq files (no id at line {})".format(i)
				print(read_id(next_line_1.strip(),next_line_2.strip()))
			elif l == 1:
				print(next_line_1.strip()+'NNN'+next_line_2.strip())
			elif l == 2:
				print('+')
			elif l == 3:
				print(next_line_1.strip()+'!!!'+next_line_2.strip())

	elif fasta:
		prev_read_1 = ''
		prev_read_2 = ''
		for next_line_1, next_line_2 in zip(reads_1, reads_2):
			if next_line_1.startswith('>'):
				print(prev_read_1+'NNN'+prev_read_2)
				prev_read_1 = ''
				prev_read_2 = ''
				print(read_id(next_line_1.strip(),next_line_2.strip()))
			else:
				prev_read_1 += next_line_1.strip()
				prev_read_2 += next_line_2.strip()

	if reads_1.readline() != '' or reads_2.readline() != '':
		print('[prophyle_paired_reads] Warning: files of different length (merged till the end of the shortest)',
				file=sys.stderr)

	reads_1.close()
	reads_2.close()

def main():
	desc="""\
		Program: prophyle_paired_end.py

		Merge paired-end FASTA or FASTQ files (possibly in gzip format) and print to stdout.
	"""
	parser = argparse.ArgumentParser(
		description=desc)

	parser.add_argument('reads_1',
						type = str,
						help = '1st FASTA or FASTQ file')
	parser.add_argument('reads_2',
						type = str,
						help = '2nd FASTA or FASTQ file')

	args = parser.parse_args()
	merge_reads(args.reads_1, args.reads_2)


if __name__ == '__main__':
	main()
