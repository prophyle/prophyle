#! /usr/bin/env python3

"""Main Prophyle file.

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT

Example:

	Download sequences:
	
		$ prophyle download bacteria
	
	Create an index for k=10 and the small testing bacterial tree:
	
		$ prophyle index -k 10 ~/prophyle/test_bacteria.nw ~/prophyle/test_viruses.nw test_idx
	
	Classify some reads:

		$ prophyle classify test_idx reads.fq > result.sam
	
Todo:
	* save configuration (trees, k, etc.) into a json; if anything changed from the last time, remove all marks
	* _is_complete should be combined with a test of files: is_missing => remove mark
	* index: automatically decide about paths for bwa, etc. (package vs. git repo)
	* index: kmer annotation to the tree
	* classificaton: support for c2, h2
	* check if prophyle_assembler & prophyle-index are newer than their source files
"""


import argparse
import datetime
import multiprocessing
import os
import shutil
import subprocess
import sys
import textwrap
import glob
import re

from . import version

c_d=os.path.dirname(os.path.realpath(__file__))
tree_d=os.path.join(c_d,"trees")

#bin_dir=os.path.dirname(__file__)
bwa=os.path.join(c_d,"prophyle-index","bwa","bwa")
ind=os.path.join(c_d,"prophyle-index","prophyle-index")
asm=os.path.join(c_d,"prophyle_assembler","prophyle_assembler")

## todo: decide about the paths for programs (execution from repo vs from package):
#    newick2makefile=os.path.join(c_d,"newick2makefile.py")
#     vs.
#    newick2makefile="prophyle_propagation_makefile.py"


newick2makefile="prophyle_propagation_makefile.py"
test_tree="prophyle_test_tree.py"
merge_fastas="prophyle_merge_fa.py"
merge_trees="prophyle_merge_trees.py"
assign="prophyle_assignment.py"

DEFAULT_K=31
DEFAULT_THREADS=multiprocessing.cpu_count()
#DEFAULT_THREADS=1
DEFAULT_MEASURE='h1'
DEFAULT_OUTPUT_FORMAT='sam'
DEFAULT_HOME_DIR=os.path.join(os.path.expanduser('~'),'prophyle')

LIBRARIES=['bacteria', 'viruses', 'plasmids', 'hmp']

FTP_NCBI='https://ftp.ncbi.nlm.nih.gov'


def _test_files(*fns,test_nonzero=False):
	"""Test if given files exist, and possibly if they are non-empty. If not, stop the program.

	Args:
		*fns: Files.
		test_nonzero (bool): Test if files have size greater than zero.

	Raises:
		AssertionError: File does not exist or it is empty.
	"""
	#print(fns)
	for fn in fns:
		assert os.path.isfile(fn), 'File "{}" does not exist'.format(fn)
		if test_nonzero:
			assert _file_sizes(fn)[0], 'File "{}" has size 0'.format(fn)


def _test_tree(fn):
	"""Test if given tree is valid for ProPhyle.

	Args:
		fn (str): Newick/NHX tree.
	"""
	_test_files(fn)
	cmd=[test_tree, '-n', fn]


def _file_sizes(*fns):
	"""Get file sizes in Bytes.

	Args:
		fns (str): File names.

	Returns:
		tuple(int): File sizes.
	"""
	return tuple( [os.stat(fn).st_size for fn in fns] )


def _run_safe(command, output_fn=None, output_fo=None):
	"""Get file sizes in Bytes.

	Args:
		output_fn (str): Name of a file for storing the output.
		output_fo (fileobject): Output file object. If both params are None, the standard output is used.

	Raises:
		RuntimeError: Command exited with non-zero code.
	"""
	assert output_fn is None or output_fo is None
	command_str=" ".join(map(lambda x: str(x),command))
	_message("Running:", command_str)
	if output_fn is None:
		if output_fo is None:
			out_fo=sys.stdout
		else:
			out_fo=output_fo
	else:
		out_fo=open(output_fn,"w+")

	error_code=subprocess.call("/bin/bash -e -o pipefail -c '{}'".format(command_str), shell=True, stdout=out_fo)

	out_fo.flush()

	if output_fn is not None:
		out_fo.close()

	if error_code==0:
		_message("Finished:", command_str)
	elif error_code==141:
		pass
		#print("Exited before finishing:", command_str, file=sys.stderr)
	else:
		_message("Unfinished, an error occurred (error code {}):".format(error_code), command_str)
		raise RuntimeError("Command error.")


def _message(*msg, upper=False):
	"""Print a ProPhyle message to stderr.

	Args:
		*msg: Message.
	"""

	dt=datetime.datetime.now()
	fdt=dt.strftime("%Y-%m-%d %H:%M:%S")
	if upper:
		msg=map(str,msg)
		msg=map(str.upper,msg)
	print('[prophyle]', fdt, " ".join(msg), file=sys.stderr)


def _touch(*fns):
	"""Touch files.

	Args:
		*fns: Files.
	"""
	for fn in fns:
		if os.path.exists(fn):
			os.utime(fn, None)
		else:
			with open(fn, 'a'):
				pass


def _rm(*fns):
	"""Remove files (might not exists).

	Args:
		*fns: Files.
	"""
	for fn in fns:
		try:
			os.remove(fn)
		except FileNotFoundError:
			pass

def _cp_to_file(fn0, fn):
	# keep rewriting attributes
	shutil.copyfile(fn0, fn)

def _cp_to_dir(fn0, d):
	# keep rewriting attributes
	shutil.copy(fn0, d)


def _makedirs(*ds):
	"""Make dirs recursively.

	Args:
		*ds: Dirs to create.
	"""
	for d in ds:
		if not os.path.isdir(d):
			cmd=['mkdir', '-p', d]
			_run_safe(cmd)


def _compile_prophyle_bin():
	"""Compile ProPhyle binaries if they don't exist yet.
	"""
	files_to_check=[
			os.path.join(c_d,'prophyle_assembler','prophyle_assembler'),
			os.path.join(c_d,'prophyle-index','prophyle-index'),
			os.path.join(c_d,'prophyle-index','bwa','bwa'),
		]
	for x in files_to_check:
		if not os.path.isfile(x):
			_message("Binaries are missing, going to compile them")
			command=["make","-C",c_d]
			_run_safe(command, output_fo=sys.stderr)
			return


def _existing_and_newer(fn0, fn):
	"""Test if file fn exists and is newer than fn0. Raise an exception if fn0 does not exist.

	Args:
		fn0 (str): Old file.
		fn (str): New file (to be generated from fn0).
	"""

	assert os.path.isfile(fn0), "Dependency '{}' does not exist".format(fn0)

	if not os.path.isfile(fn):
		return False

	if os.path.getmtime(fn0)<=os.path.getmtime(fn):
		return True
	else:
		return False


#####################
# PROPHYLE DOWNLOAD #
#####################

def __mark_fn(d, i, name):
	"""Create a mark name.

	Args:
		d (str): Directory.
		i (int): Number of the step.
		name (str): Name of the mark.
	"""
	if name is None:
		return os.path.join(d,".complete.{}".format(i))
	else:
		return os.path.join(d,".complete.{}.{}".format(name,i))


def _mark_complete(d, i=1, name=None):
	"""Create a mark file (an empty file to mark a finished step nb i).

	Args:
		d (str): Directory.
		i (int): Number of the step.
		name (str): Name of the mark.
	"""

	assert i>0

	_touch(__mark_fn(d, i, name))

 
def _is_complete(d, i=1, name=None):
	"""Check if a mark file i exists AND is newer than the mark file (i-1).

	Args:
		d (str): Directory.
		i (int): Number of the step.
		name (str): Name of the mark.
	"""

	assert i>0
	fn=__mark_fn(d,i, name)
	fn0=__mark_fn(d,i-1, name)

	if i==1:
		return os.path.isfile(fn)
	else:
		return _existing_and_newer(fn0, fn)


def _missing_library(d):
	"""Check if a mark file i exists AND is newer than the mark file (i-1). Create the library dir.

	Args:
		d (str): Directory.
		i (int): Number of the step.
		name (str): Name of the mark.
	"""

	l=os.path.dirname(d)
	_makedirs(d)
	if _is_complete(d,1):
		_message("Skipping downloading library '{}' (already exists)".format(l))
		return False
	else:
		_message("Downloading library '{}'".format(l))
		return True


def _pseudo_fai(d):
	"""Generate a psedudofai file for given directory (directory/*.fa => directory.fai).

	Pseudofai format = TSV with 2 two columns: filename, sequence header (text after > in FASTA).

	Args:
		d (str): Directory.
	"""
	l=os.path.dirname(d)
	pseudofai_fn=d+".pseudofai"
	_makedirs(d)
	if _is_complete(d,2) and os.path.isfile(pseudofai_fn):
		_message("Skipping generating pseudofai for library '{}' (already exists)".format(l))
	else:
		_message("Generating pseudofai for library '{}'".format(l))
		assert d[-1]!="/"
		#cmd=['grep -r --include=\\*.{fa,ffn,fna}', '">"', d, '| sed "s/:>/\t/"']
		cmd=[
			'find', d, '-name', "'*.fa'", "-o", "-name", "'*.ffn'", "-o", "-name", "'*.fna'", "-exec", "grep", "-H", '">"', "{}", "\\;",
			"|", 'sed', '"s/\:>/\t/"']

		_run_safe(cmd, pseudofai_fn)
		_mark_complete(d, 2)


def prophyle_download(library, library_dir):
	"""Create a library Download genomic library and copy the corresponding tree.

	Args:
		library (str): Library to download (bacteria / viruses / ...)
		library_dir (str): Directory where download files will be downloaded.

	Todo:
		* Add support for alternative URLs (http / ftp, backup refseq sites, etc.).
			* http://downloads.hmpdacc.org/data/HMREFG/all_seqs.fa.bz2
			* ftp://public-ftp.hmpdacc.org/HMREFG/all_seqs.fa.bz2
	"""

	if library=="all":
		for l in LIBRARIES:
			download(l, library_dir)
		return
	else:
		assert library in LIBRARIES

	if library_dir is None:
		d=os.path.join(os.path.expanduser("~/prophyle"),library)
	else:
		d=os.path.join(library_dir,library)
	#print('making',d, file=sys.stderr)
	#os.makedirs(d, exist_ok=True)
	_makedirs(d)

	_message("Checking library '{}' in '{}'".format(library,d))

	lib_missing=_missing_library(d)
	if lib_missing:
		for test_prefix in ["","test_"]:
			fn="{}{}.nw".format(test_prefix,library,)
			nhx=os.path.join(tree_d,fn)
			new_nhx=os.path.join(d,"..",fn)
			_test_files(nhx)
			_message("Copying Newick/NHX tree '{}' to '{}'".format(nhx,new_nhx))
			_cp_to_file(nhx, new_nhx)

	if library=='bacteria':
		if lib_missing:
			cmd=['cd', d, '&& curl', FTP_NCBI+'/genomes/archive/old_refseq/Bacteria/all.fna.tar.gz | tar xz']
			_run_safe(cmd)
			_mark_complete(d, 1)
		#_pseudo_fai(d)

	elif library=='viruses':
		if lib_missing:
			#cmd=['cd', d, '&& curl', FTP_NCBI+'/genomes/Viruses/all.ffn.tar.gz | tar xz']
			#_run_safe(cmd)
			cmd=['cd', d, '&& curl', FTP_NCBI+'/genomes/Viruses/all.fna.tar.gz | tar xz']
			_run_safe(cmd)
			_mark_complete(d, 1)
		#_pseudo_fai(d)

	elif library=='plasmids':
		if lib_missing:
			cmd=['cd', d, '&& curl', FTP_NCBI+'/genomes/archive/old_refseq/Plasmids/plasmids.all.fna.tar.gz | tar xz --strip 5']
			_run_safe(cmd)
			_mark_complete(d, 1)
		#_pseudo_fai(d)

	elif library=='hmp':
		if lib_missing:
			# fix when error appears
			cmd=['cd', d, '&& curl http://downloads.hmpdacc.org/data/HMREFG/all_seqs.fa.bz2 | bzip2 -d']
			_run_safe(cmd,os.path.join(d,"all_seqs.fa"))
			_mark_complete(d, 1)
		#_pseudo_fai(d)

	else:
		raise ValueError('Unknown library ""'.format(library))


##################
# PROPHYLE INDEX #
##################

def _create_makefile(index_dir, k, library_dir):
	"""Create a Makefile for k-mer propagation.

	Args:
		index_dir (str): Index directory.
		k (int): K-mer size.
		library_dir (library_dir): Library directory.

	Todo:
		* Add checking of params.mk
	"""
	_message('Creating Makefile for k-mer propagation')
	propagation_dir=os.path.join(index_dir, 'propagation')
	_makedirs(propagation_dir)

	makefile=os.path.join(propagation_dir,'Makefile')
	tree_fn=os.path.join(index_dir,'tree.nw')
	_test_tree(tree_fn)
	#_test_files(newick2makefile, tree_fn)
	command=[newick2makefile, '-n', tree_fn, '-k', k, '-o', './', '-l', os.path.abspath(library_dir)]

	with open(os.path.join(propagation_dir, "params.mk"),"w+") as f:
		f.write("K={}\n".format(k))
	_run_safe(command,makefile)


def _propagate(index_dir,threads):
	"""Run k-mer propagation.

	Args:
		index_dir (str): Index directory.
		threads (int): Number of threads for Makefile.
	"""
	_message('Running k-mer propagation')
	propagation_dir=os.path.join(index_dir, 'propagation')
	_test_files(os.path.join(propagation_dir, 'Makefile'),test_nonzero=True)
	command=['make', '-j', threads, '-C', propagation_dir, 'V=1', "PRG_ASM={}".format(asm)]
	_run_safe(command)


def _merge_trees(in_trees, out_tree, no_prefixes):
	"""Merge input trees into a single tree.

	Args:
		in_trees (list of str): Input NHX trees.
		out_tree (str): Output NHX tree.
	"""

	_message('Generating index tree')
	_test_files(*in_trees)
	command=[merge_trees] + in_trees + [out_tree]
	if no_prefixes:
		command += ['-P']
	_run_safe(command)


def _merge_fastas(index_dir):
	"""Merge reduced FASTA files after k-mer propagation and create index.fa.

	Args:
		index_dir (str): Index directory.

	Todo:
		* check files for all nodes exist and are of size > 0
	"""

	_message('Generating index.fa')
	propagation_dir=os.path.join(index_dir, 'propagation')
	index_fa=os.path.join(index_dir,"index.fa")
	#_test_files(merge_fastas)
	command=[merge_fastas, propagation_dir]
	_run_safe(command, index_fa)
	_touch(index_fa+".complete")


def _fa2pac(fa_fn):
	"""Run `bwa fa2pac` (FA => 2bit).

	Args:
		fa_fn (str): FASTA file.
	"""

	_message('Generating packed FASTA file')
	_test_files(bwa, fa_fn)
	command=[bwa, 'fa2pac', fa_fn, fa_fn]
	_run_safe(command)


def _pac2bwt(fa_fn):
	"""Run `bwa pac2bwtgen` (2bit => BWT).

	Args:
		fa_fn (str): FASTA file.
	"""

	_message('Generating BWT')
	_test_files(bwa, fa_fn+".pac")
	command=[bwa, 'pac2bwtgen', fa_fn+".pac", fa_fn+".bwt"]
	_run_safe(command)


def _bwt2bwtocc(fa_fn):
	"""Run `bwa bwtupdate` (BWT => BWT+OCC).

	Args:
		fa_fn (str): FASTA file.
	"""

	_message('Generating sampled OCC array')
	_test_files(bwa, fa_fn+".bwt")
	command=[bwa, 'bwtupdate', fa_fn+".bwt"]
	_run_safe(command)


def _bwtocc2sa(fa_fn):
	"""Run `bwa bwt2sa` (BWT+OCC => SSA).

	Args:
		fa_fn (str): FASTA file.
	"""

	_message('Generating sampled SA')
	_test_files(bwa, fa_fn+".bwt")
	command=[bwa, 'bwt2sa', fa_fn+".bwt", fa_fn+".sa"]
	_run_safe(command)


def _bwtocc2klcp(fa_fn,k):
	"""Create k-LCP `` (BWT => k-LCP).

	Args:
		fa_fn (str): FASTA file.
		k (int): K-mer size.
	"""

	_message('Generating k-LCP array')
	_test_files(ind, fa_fn+".bwt")
	command=[ind, 'build', '-k', k, fa_fn]
	_run_safe(command)


def _bwtocc2sa_klcp(fa_fn,k):
	"""Create k-LCP `` (BWT => k-LCP).

	Args:
		fa_fn (str): FASTA file.
		k (int): K-mer size.
	"""

	_message('Generating k-LCP array and SA in parallel')
	_test_files(ind, fa_fn+".bwt")
	command=[ind, 'build', '-s', '-k', k, fa_fn]
	_run_safe(command)


def prophyle_index(index_dir, threads, k, trees_fn, library_dir, construct_klcp, force, no_prefixes):
	"""Build a Prophyle index.

	Args:
		index_dir (str): Index directory.
		threads (int): Number of threads in k-mer propagation.
		k (int): K-mer size.
		tree_fn (str): Newick/NHX tree.
		library_dir (str): Library directory.
		klcp (bool): Generate klcp.
		force (bool): Rewrite files if they already exist.
		no_prefixes (bool): Don't prepend prefixes to node names during tree merging.

	Todo:
		* klcp in parallel with SA
		* copy Newick only if it is newer
		* add update the tree with number of k-mers
	"""

	assert isinstance(k, int)
	assert isinstance(threads, int)
	assert k>1
	assert threads>0

	_compile_prophyle_bin()


	index_fa=os.path.join(index_dir,'index.fa')
	index_tree=os.path.join(index_dir,'tree.nw')
	makefile_dir=os.path.join(index_dir,'propagation')
	makefile=os.path.join(index_dir,'propagation','Makefile')

	# recompute = recompute everything from now on
	# force==True => start to recompute everything from beginning
	recompute=force

	# make index dir
	_makedirs(index_dir)

	#
	# 1) Newick
	#
	
	#if not _existing_and_newer(tree_fn, index_tree):
	if not _is_complete(index_dir, 1):
		recompute=True


	if recompute:
		_message('[1/5] Copying/merging trees', upper=True)
		_merge_trees(trees_fn, index_tree, no_prefixes=no_prefixes)
		_mark_complete(index_dir, 1)
	else:
		_message('[1/5] Tree already exists, skipping copying', upper=True)

	#
	# 2) Create and run Makefile for propagation, and merge FASTA files
	#

	if not _is_complete(index_dir, 2):
		recompute=True

	if recompute:
		# todo: check if something should be deleted (e.g., the propagation dir)
		_message('[2/5] Running k-mer propagation', upper=True)
		_create_makefile(index_dir, k, library_dir)
		_propagate(index_dir, threads=threads)
		_merge_fastas(index_dir)
		_mark_complete(index_dir, 2)
	else:
		_message('[2/5] K-mers have already been propagating, skipping propagation', upper=True)

	#
	# 3) BWT + OCC
	#

	if not _is_complete(index_dir, 3):
		recompute=True

	#if ccontinue and os.path.isfile(index_fa+'.bwt') and os.path.isfile(index_fa+'.bwt.complete'):

	if recompute:
		_message('[3/5] Constructing BWT+OCC', upper=True)
		_rm(index_fa+'.bwt',index_fa+'.bwt.complete')
		_fa2pac(index_fa)
		_pac2bwt(index_fa)
		_bwt2bwtocc(index_fa)
		_mark_complete(index_dir, 3)
	else:
		_message('[3/5] BWT and OCC already exist, skipping their construction', upper=True)

	#
	# 4) SA + 5) KLCP (compute SA + KLCP in parallel)
	#

	klcp_fn="{}.{}.klcp".format(index_fa,k)

	if construct_klcp:

		if not _is_complete(index_dir, 4):
			#SA not computed yet => compute it in parallel with KLCP
			recompute=True

		if recompute:
			_message('[4/5],[5/5] Constructing SA + KLCP in parallel ', upper=True)
			_bwtocc2sa_klcp(index_fa, k)
			_mark_complete(index_dir, 4)
			_mark_complete(index_dir, 5)
			return

	#
	# 4) SA (compute only SA)
	#

	if not _is_complete(index_dir, 4):
		recompute=True

	if recompute:
		_message('[4/5] Constructing SA', upper=True)
		_bwtocc2sa(index_fa)
	else:
		_message('[4/5] SA already exists, skipping its construction', upper=True)


	#
	# 5) KLCP (compute only KLCP)
	#

	if construct_klcp:
		if not _is_complete(index_dir, 5):
			recompute=True

		if recompute:
			_message('[5/5] Constructing k-LCP', upper=True)
			_bwtocc2klcp(index_fa,k)
			_mark_complete(index_dir, 5)
		else:
			_message('[5/5] k-LCP already exists, skipping its construction', upper=True)


#####################
# PROPHYLE CLASSIFY #
#####################

def prophyle_classify(index_dir,fq_fn,k,use_rolling_window,out_format,mimic_kraken,measure,annotate,tie_lca):

	"""Run Prophyle classification.

	Args:
		index_dir (str): Index directory.
		fq_fn (str): Input reads.
		k (int): K-mer size (None => detect automatically).
		use_rolling_window (bool): Use rolling window.
		out_format (str): Output format: sam / kraken.
		mimic_kraken (bool): Mimic Kraken algorithm (compute LCA for each k-mer).
		measure (str): Measure used for classification (h1 / h2 / c1 / c2).
		annotate (bool): Annotate assignments (insert annotations from Newick to SAM).
		tie_lca (bool): If multiple equally good assignments found, compute their LCA.
	"""


	_compile_prophyle_bin()
	index_fa=os.path.join(index_dir, 'index.fa')
	index_tree=os.path.join(index_dir, 'tree.nw')

	if k is None:
		klcps=glob.glob(os.path.join(index_dir,"*.klcp"))

		assert len(klcps)<2, "K-mer length could not be detected (several k-LCP files exist). Please use the '-k' parameter."
		assert len(klcps)>0, "K-mer length could not be detected (no k-LCP file exists). Please use the '-k' parameter."
		klcp=klcps[0]

		re_klcp=re.compile(r'.*/index\.fa\.([0-9]+)\.klcp$')
		klcp_match=re_klcp.match(klcp)
		k=klcp_match.group(1)
		_message("Automatic detection of k-mer length: k={}".format(k))

	_test_tree(index_tree)
	#_test_files(fq_fn,index_fa,ind,assign)

	if fq_fn!="-":
		_test_files(fq_fn)

	_test_files(index_fa,ind)

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

	if use_rolling_window:
		klcp_fn="{}.{}.klcp".format(index_fa,k)
		_test_files(klcp_fn)
		(klcp_s,)=_file_sizes(klcp_fn)
		assert abs(bwt_s - 4*klcp_s) < 1000, 'Inconsistent index (KLCP vs. BWT)'

	if mimic_kraken:
		cmd_assign=[assign, '-i', '-', '-k', k, '-n', index_tree, '-m', 'h1', '-f', 'kraken', '-l', '-t']
	else:
		cmd_assign=[assign, '-i', '-', '-k', k, '-n', index_tree, '-m', measure, '-f', out_format]
		if annotate:
			cmd_assign+=['--annotate']
		if tie_lca:
			cmd_assign+=['--tie-lca']

	cmd_query=[ind, 'query', '-k', k, '-u' if use_rolling_window else '', index_fa, fq_fn]


	#(['|', '|'] if mimic_kraken else ['|']) \
	command=cmd_query + ['|'] + cmd_assign
	_run_safe(command)


########
# MAIN #
########

def parser():

	class MyParser(argparse.ArgumentParser):
		def error(self, message):
			if len(sys.argv)==2:
				self.print_help()
			else:
				print('error: {}'.format(message), file=sys.stderr)
			sys.exit(2)

	desc="""\
		Program: prophyle (phylogeny-based metagenomic classification)
		Version: {V}
		Authors: Karel Brinda <kbrinda@hsph.harvard.edu>, Kamil Salikhov <kamil.salikhov@univ-mlv.fr>,
		         Simone Pignotti <pignottisimone@gmail.com>, Gregory Kucherov <gregory.kucherov@univ-mlv.fr>

		Usage:   prophyle <command> [options]
	""".format(V=version.VERSION)
	parser = MyParser(formatter_class=argparse.RawDescriptionHelpFormatter,description=textwrap.dedent(desc))
	subparsers = parser.add_subparsers(help="",description=argparse.SUPPRESS,dest='subcommand',metavar="")
	fc=lambda prog: argparse.HelpFormatter(prog,max_help_position=27)

	##########

	parser_download = subparsers.add_parser('download',
			help='download a genomic database',
			#description='Download RefSeq and HMP databases.',
			formatter_class=fc,
		)
	parser_download.add_argument(
			'library',
			metavar='<library>',
			choices=LIBRARIES+['all'],
			help='genomic library {}'.format(LIBRARIES+['all']),
		)
	parser_download.add_argument(
			'-d',
			metavar='DIR',
			dest='home_dir',
			type=str,
			default=None,
			help='directory for the tree and the sequences [~/prophyle]',
		)

	##########

	parser_index = subparsers.add_parser('index',
			help='build index',
			#description='Build a ProPhyle index (i.e., propagate k-mers and construct a BWT-index with k-LCP).',
			formatter_class=fc,
		)
	parser_index.add_argument('tree',
			metavar='<tree.nw>',
			type=str,
			nargs='+',
			help='phylogenetic tree (in Newick/NHX)',
		)
	parser_index.add_argument(
			'index_dir',
			metavar='<index.dir>',
			type=str,
			help='index directory (will be created)',
		)
	parser_index.add_argument(
			'-g',
			metavar='DIR',
			dest='library_dir',
			type=str,
			help='directory with the library sequences [directory of the first tree]',
			default=None,
			#required=True,
		)
	parser_index.add_argument(
			'-j',
			metavar='INT',
			dest='threads',
			type=int,
			help='number of threads [auto ({})]'.format(DEFAULT_THREADS),
			default=DEFAULT_THREADS,
		)
	parser_index.add_argument(
			'-k',
			dest='k',
			metavar='INT',
			type=int,
			help='k-mer length [{}]'.format(DEFAULT_K),
			default=DEFAULT_K,
		)
	parser_index.add_argument(
			'-F',
			dest='force',
			action='store_true',
			help='rewrite index files if they already exist',
		)
	parser_index.add_argument(
			'-P',
			dest='no_prefixes',
			action='store_true',
			help='do not add prefixes to node names when multiple trees are used',
		)
	parser_index.add_argument(
			'-K',
			dest='klcp',
			action='store_false',
			help='skip k-LCP construction',
		)

	##########

	parser_classify = subparsers.add_parser('classify',
			help='classify reads',
			#description='Classify reads.',
			formatter_class=fc,
		)
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
			help='file with reads in FASTA or FASTQ (use - for standard input)',
		)
	parser_classify.add_argument(
			'-k',
			dest='k',
			metavar='INT',
			type=int,
			help='k-mer length [detect automatically]',
			default=None,
		)
	parser_classify.add_argument(
			'-R',
			dest='rolling_window',
			action='store_false',
			help='use restarted search for matching rather than rolling window (slower, but k-LCP is not needed)',
		)
	parser_classify.add_argument(
			'-m',
			dest='measure',
			choices=['h1','c1'],
			help='measure: h1=hit count, c1=coverage [{}]'.format(DEFAULT_MEASURE),
			default=DEFAULT_MEASURE,
		)
	parser_classify.add_argument(
			'-f',
			dest='oform',
			choices=['kraken','sam'],
			default=DEFAULT_OUTPUT_FORMAT,
			help='output format [{}]'.format(DEFAULT_OUTPUT_FORMAT),
		)
	parser_classify.add_argument(
			'-A',
			dest='annotate',
			action='store_true',
			help='annotate assignments',
		)
	parser_classify.add_argument(
			'-L',
			dest='tie',
			action='store_true',
			help='use LCA when tie (multiple hits with the same score)',
		)
	parser_classify.add_argument(
			'-M',
			dest='mimic',
			action='store_true',
			help='mimic Kraken algorithm and output (for debugging purposes)',
		)

	##########

	return parser

def main():
	try:
		par=parser()
		args = par.parse_args()
		subcommand=args.subcommand

		if subcommand=="download":
			prophyle_download(
					library=args.library,
					library_dir=args.home_dir,
				)
			_message('Downloading finished')

		elif subcommand=="index":
			if args.library_dir is None:
				library_dir=os.path.dirname(args.tree[0])
			else:
				library_dir=args.library_dir
			prophyle_index(
					index_dir=args.index_dir,
					threads=args.threads,
					k=args.k,
					trees_fn=args.tree,
					library_dir=library_dir,
					force=args.force,
					construct_klcp=args.klcp,
					no_prefixes=args.no_prefixes,
				)
			_message('Index construction finished')

		elif subcommand=="classify":
			prophyle_classify(
					index_dir=args.index_dir,
					fq_fn=args.reads,
					k=args.k,
					use_rolling_window=args.rolling_window,
					out_format=args.oform,
					mimic_kraken=args.mimic,
					measure=args.measure,
					tie_lca=args.tie,
					annotate=args.annotate,
				)
			_message('Classificaton finished')

		else:
			msg_lns=par.format_help().split("\n")[2:]
			msg_lns=[x for x in msg_lns if x.find("optional arguments")==-1 and x.find("show this help")==-1]
			msg="\n".join(msg_lns)
			msg=msg.replace("\n\n",'\n').replace("subcommands:\n","Command:").replace("Usage","\nUsage")
			print(file=sys.stderr)
			print(msg,file=sys.stderr)
			sys.exit(1)

	except BrokenPipeError:
		# pipe error (e.g., when head is used)
		sys.stderr.close()
		exit(0)

if __name__ == "__main__":
	main()
