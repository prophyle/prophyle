#! /usr/bin/env python3

import argparse
import datetime
import multiprocessing
import os
import shutil
import subprocess
import sys

#bin_dir=os.path.normpath(os.path.join(os.path.dirname(__file__),"../../bin"))
bin_dir=os.path.dirname(__file__)
bwa=os.path.join(bin_dir,"bwa")
exk=os.path.join(bin_dir,"exk")
asm=os.path.abspath(os.path.join(bin_dir,"assembler"))
newick2makefile=os.path.join(bin_dir,"newick2makefile.py")
test_newick=os.path.join(bin_dir,"test_newick_tree.py")
merge_fastas=os.path.join(bin_dir,"create_final_fasta.py")
assign=os.path.join(bin_dir,"assignment.py")

DEFAULT_K=32
DEFAULT_THREADS=multiprocessing.cpu_count()
DEFAULT_MEASURE='h1'

def _test_files(*fns):
	#print(fns)
	for fn in fns:
		assert os.path.isfile(fn), 'File "{}" does not exist'.format(fn)

def _file_sizes(*fns):
	return (os.stat(fn).st_size for fn in fns)

def _run_safe(command, output_fn=None):
	command_str=" ".join(map(lambda x: str(x),command))
	print("Running:", command_str, file=sys.stderr)
	if output_fn is None:
		out_fo=sys.stdout
	else:
		out_fo=open(output_fn,"w+")
	error_code=subprocess.call("/bin/bash -x -o pipefail -c '{}'".format(command_str), shell=True, stdout=out_fo)
	if error_code==0:
		print("Finished:", command_str, file=sys.stderr)
	else:
		print("Finished with error (error code {}):".format(error_code), command_str, file=sys.stderr)

def _message(msg):
	dt=datetime.datetime.now()
	fdt=dt.strftime("%Y-%m-%d %H:%M:%S")
	print('[metang]', fdt, msg, file=sys.stderr)

###############
# METANG INIT #
###############

def init():
	pass


################
# METANG INDEX #
################

def _create_makefile(index_dir, k, library_dir):
	_message('Creating Makefile for k-mer propagation')
	propagation_dir=os.path.join(index_dir, 'propagation')
	os.makedirs(propagation_dir)

	makefile=os.path.join(propagation_dir,'Makefile')
	newick_fn=os.path.join(index_dir,'tree.newick')
	#_test_files(newick2makefile, newick_fn)
	command=[newick2makefile, '-n', newick_fn, '-k', k, '-o', './', '-l', os.path.abspath(library_dir)]
	_run_safe(command,makefile)

def _propagate(index_dir,threads):
	_message('Running k-mer propagation')
	propagation_dir=os.path.join(index_dir, 'propagation')
	_test_files(os.path.join(propagation_dir, 'Makefile'))
	command=['make', '-j', threads, '-C', propagation_dir, 'V=1', "ASSEMBLER={}".format(asm)]
	_run_safe(command)

def _merge_fastas(index_dir):
	_message('Generating index.fa')
	propagation_dir=os.path.join(index_dir, 'propagation')
	# todo: check files for all nodes exist and are of size > 0
	index_fa=os.path.join(index_dir,"index.fa")
	_test_files(merge_fastas)
	command=[merge_fastas, propagation_dir]
	_run_safe(command, index_fa)	

def _fa2pac(fa_fn):
	_message('Generating packed FASTA file')
	_test_files(bwa, fa_fn)
	command=[bwa, 'fa2pac', fa_fn, fa_fn]
	_run_safe(command)

def _pac2bwt(fa_fn):
	_message('Generating BWT')
	_test_files(bwa, fa_fn+".pac")
	command=[bwa, 'pac2bwtgen', fa_fn+".pac", fa_fn+".bwt"]
	_run_safe(command)

def _bwt2bwtocc(fa_fn):
	_message('Generating sampled OCC array')
	_test_files(bwa, fa_fn+".bwt")
	command=[bwa, 'bwtupdate', fa_fn+".bwt"]
	_run_safe(command)

def _bwtocc2sa(fa_fn):
	_message('Generating sampled SA')
	_test_files(bwa, fa_fn+".bwt")
	command=[bwa, 'bwt2sa', fa_fn+".bwt", fa_fn+".sa"]
	_run_safe(command)

def _bwtocc2klcp(fa_fn,k):
	_message('Generating k-LCP array')
	_test_files(exk, fa_fn+".bwt")
	command=[exk, 'index', '-k', k, fa_fn]
	_run_safe(command)

def index(index_dir, threads, k, newick_fn, library_dir, cont=False, klcp=True):
	assert k>1

	# check files & dirs
	_test_files(newick_fn)
	index_fa=os.path.join(index_dir,'index.fa')
	index_newick=os.path.join(index_dir,'tree.newick')
	makefile_dir=os.path.join(index_dir,'propagation')
	makefile=os.path.join(index_dir,'propagation','Makefile')

	assert not os.path.isfile(index_dir) 
	assert not os.path.isdir(index_dir)

	# make index dir
	os.makedirs(index_dir)

	# copy newick
	shutil.copy(newick_fn, index_newick)

	# create makefile
	_create_makefile(index_dir, k, library_dir)

	# run makefile
	_propagate(index_dir, threads=threads)

	# merge fastas
	_merge_fastas(index_dir)

	# bwa index & klcp
	_fa2pac(index_fa)
	_pac2bwt(index_fa)
	_bwt2bwtocc(index_fa)
	_bwtocc2sa(index_fa)
	if klcp:
		_bwtocc2klcp(index_fa,k)


###################
# METANG CLASSIFY #
###################

def classify(index_dir,fq_fn,k,use_klcp,out_format,mimic_kraken,measure,annotate,tie_lca):
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

	(bwt_s, sa_s, pac_s)=_file_sizes(index_fa+'.bwt',index_fa+'.sa',index_fa+'.pac')
	assert abs(bwt_s - 2*sa_s) < 1000, 'Inconsistent index (SA vs. BWT)'
	assert abs(bwt_s - 2*pac_s) < 1000, 'Inconsistent index (PAC vs. BWT)'

	if use_klcp:
		klcp_fn="{}.{}.bit.klcp".format(index_fa,k)
		_test_files(klcp_fn)
		(klcp_s,)=_file_sizes(klcp_fn)
		assert abs(bwt_s - 2*klcp_s) < 1000, 'Inconsistent index (KLCP vs. BWT)'

	# todo: add integrity checking (correct file size: |sa|=|pac|, |bwt|=2|sa|)

	if mimic_kraken:
		cmd_assign=[assign, '-i', '-', '-k', k, '-n', index_newick, '-m', 'h1', '-f', 'kraken', '-l', '-t']
	else:
		cmd_assign=[assign, '-i', '-', '-k', k, '-n', index_newick, '-m', measure, '-f', out_format]
		if annotate:
			cmd_assign+=['--annotate']
		if tie_lca:
			cmd_assign+=['--tie-lca']

	cmd_match=[exk, 'match', '-k', k, '-u' if use_klcp else '', index_fa, fq_fn]


	#(['|', '|'] if mimic_kraken else ['|']) \
	command=cmd_match + ['|'] + cmd_assign
	_run_safe(command)


########
# MAIN #
########

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	subparsers = parser.add_subparsers(help='sub-command help',dest='subcommand')
	fc=lambda prog: argparse.HelpFormatter(prog,max_help_position=27)

	##########

	parser_init = subparsers.add_parser('init', help='Initialize data', formatter_class=fc)
	parser_init.add_argument('--bar', type=int, help='bar help', required=False)

	##########

	parser_index = subparsers.add_parser('index', help='Create index', formatter_class=fc)
	parser_index.add_argument(
			'-n','--newick',
			metavar='FILE',
			dest='newick',
			type=str,
			help='taxonomy tree (in Newick format)',
			required=True,
		)
	parser_index.add_argument(
			'index_dir',
			metavar='<index.dir>',
			type=str,
			help='index directory (will be created)',
		)
	parser_index.add_argument(
			'-l','--lib-dir',
			metavar='DIR',
			dest='library_dir',
			type=str,
			help='directory with genomic sequences',
			required=True,
		)
	parser_index.add_argument(
			'-t','--threads',
			metavar='INT',
			dest='threads',
			type=int,
			help='number of threads [auto={}]'.format(DEFAULT_THREADS),
			default=DEFAULT_THREADS,
		)
	parser_index.add_argument(
			'-k','--kmer-len',
			dest='k',
			metavar='INT',
			type=int,
			help='k-mer length [{}]'.format(DEFAULT_K),
			default=DEFAULT_K,
		)

	##########

	parser_classify = subparsers.add_parser('classify', help='Classify reads', formatter_class=fc)
	parser_classify.add_argument(
			'index_dir',
			metavar='<index.dir>',
			type=str,
			help='index directory',
		)
	parser_classify.add_argument(
			'reads',
			metavar='<reads.fq>',
			type=str,
			help='file with reads in FASTA or FASTQ [- for standard input]',
		)
	parser_classify.add_argument(
			'-k','--kmer-len',
			dest='k',
			metavar='INT',
			type=int,
			help='k-mer length [{}]'.format(DEFAULT_K),
			default=DEFAULT_K,
		)
	parser_classify.add_argument(
			'-n','--no-klcp',
			dest='klcp',
			action='store_false',
			help='do not use k-LCP',
		)
	parser_classify.add_argument(
			'-m','--measure',
			dest='measure',
			choices=['h1','c1'],
			help='measure: h1=hit count, c1=coverage [{}]'.format(DEFAULT_MEASURE),
			default=DEFAULT_MEASURE,
		)
	parser_classify.add_argument(
			'-o','--out-form',
			dest='oform',
			choices=['kraken','sam'],
			default='sam',
			help='output format',
		)
	parser_classify.add_argument(
			'--annotate',
			dest='annotate',
			action='store_true',
			help='annotate assignments',
		)
	parser_classify.add_argument(
			'--tie-lca',
			dest='tie',
			action='store_true',
			help='use LCA when tie (multiple hits with the same score)',
		)
	parser_classify.add_argument(
			'--mimic-kraken',
			dest='mimic',
			action='store_true',
			help='mimic Kraken algorithm and output (for debugging purposes)',
		)

	##########

	args = parser.parse_args()
	subcommand=args.subcommand

	if subcommand=="init":
		init()

	elif subcommand=="index":		
		index(
				index_dir=args.index_dir,
				threads=args.threads,
				k=args.k,
				newick_fn=args.newick,
				library_dir=args.library_dir,
			)

	elif subcommand=="classify":
		classify(
				index_dir=args.index_dir,
				fq_fn=args.reads,
				k=args.k,
				use_klcp=args.klcp,
				out_format=args.oform,
				mimic_kraken=args.mimic,
				measure=args.measure,
				tie_lca=args.tie,
				annotate=args.annotate,
			)

	else:
		parser.print_help()
		sys.exit(1)
