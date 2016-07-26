#! /usr/bin/env python3

import argparse
import os
import subprocess
import sys

#bin_dir=os.path.normpath(os.path.join(os.path.dirname(__file__),"../../bin"))
bin_dir=os.path.dirname(__file__)
#print(bin_dir)
bwa=os.path.join(bin_dir,"bwa")
exk=os.path.join(bin_dir,"exk")
assign=os.path.join(bin_dir,"assignment.py")

def _test_files(*fns):
	#print(fns)
	for fn in fns:
		assert os.path.isfile(fn), 'File "{}" does not exist'.format(fn)

def init():
	pass

def index():
	pass

def classify(index_dir,fq_fn,k,use_klcp,out_format='sam'):
	index_fa=os.path.join(index_dir, 'index.fa')
	index_newick=os.path.join(index_dir, 'tree.newick')

	_test_files(fq_fn,index_fa,index_newick,exk,assign)

	_test_files(
			index_fa+'.bwt',
			index_fa+'.pac',
			index_fa+'.sa',
			index_fa+'.ann',
			index_fa+'.amb',
		)

	if use_klcp:
		_test_files("{}.{}.bit.klcp".format(index_fa,k))

	# todo: add integrity checking (correct file size: |sa|=|pac|, |bwt|=2|sa|)

	command=[exk, 'match', '-k', str(k), '-u' if use_klcp else '', index_fa, fq_fn] \
		+ \
		['|'] \
		+ \
		[assign, '-i', '-', '-k', str(k), '-n', index_newick, '-f', out_format]
	subprocess.run("/bin/bash -x -o pipefail -c '{}'".format(" ".join(command)), shell=True, check=True)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()

	subparsers = parser.add_subparsers(help='sub-command help',dest='subcommand')

	parser_init = subparsers.add_parser('init', help='Initialize data')
	parser_init.add_argument('--bar', type=int, help='bar help', required=False)

	parser_index = subparsers.add_parser('index', help='Create index')
	parser_index.add_argument('--bar', type=int, help='bar help', required=False)

	parser_classify = subparsers.add_parser('classify', help='Classify reads')
	parser_classify.add_argument(
			'-i',
			metavar='index_dir',
			dest='index_dir',
			type=str,
			help='Directory with index',
			required=True,
		)
	parser_classify.add_argument(
			'-f',
			metavar='reads.fq',
			dest='reads',
			type=str,
			help='File with reads in FASTA or FASTQ [- for standard input]',
			required=True,
		)
	parser_classify.add_argument(
			'-k',
			dest='k',
			type=int,
			help='K-mer length',
			required=True,
		)
	parser_classify.add_argument(
			'--no-klcp','-n',
			dest='klcp',
			action='store_false',
			help='Do not use k-LCP',
		)

	args = parser.parse_args()
	subcommand=args.subcommand

	if subcommand=="init":
		init()
	elif subcommand=="index":

		index()
	elif subcommand=="classify":
		classify(
				index_dir=args.index_dir,
				fq_fn=args.reads,
				k=args.k,
				use_klcp=args.klcp,
			)
	else:
		parser.print_help()
		sys.exit(1)
