#! /usr/bin/env python3

import os
import shutil
import datetime
import sys
import argparse

from ete3 import Tree

import logging

DEFAULT_FORMAT = 1

#logger = logging.getLogger()
#handler = logging.StreamHandler()
#formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s')
#handler.setFormatter(formatter)
#logger.addHandler(handler)
#logger.setLevel(logging.INFO)

def size_in_mb(file_fn):
	return os.path.getsize(file_fn)/(1024**2)

def merge_fasta_files(input_files,output_file,is_leaf):
	"""Merge files, remove empty lines.
	"""

#	if len(input_files)==1:
#
#		ln_name=os.path.basename(output_file)
#		ln_dir=os.path.dirname(output_file)
#		rel_path=os.path.relpath(input_files[0],ln_dir)
#
#		cmd =  (
#				"{o}: {i}\n" +
#				#"{o}:\n" +
#				"\t(cd {d} && ln -sf {i2} {o2})\n" +
#				#"\ttouch {o}\n"
#				"\n"
#			).format(
#				o=output_file,
#				i=input_files[0],
#				d=ln_dir,
#				i2=rel_path,
#				o2=ln_name,
#			)
#
#
#	else:
	if is_leaf:
		cmd =  (
				"{o}: {i}\n" +
				"\tcat $^ | $(MASKING) > $@\n\n"
			).format(
				i=' \\\n\t\t'.join(input_files),
				o=output_file,
			)
	else:
		cmd =  (
					"{o}: {i}\n" +
					"\tcat $^ > $@\n\n"
				).format(
					i=' \\\n\t\t'.join(input_files),
					o=output_file,
				)
	print(cmd)

def assembly(input_files, output_files, intersection_file):
	assert(len(input_files)==len(output_files))
	#logger.info('Starting assembly. Input files: {}. Output files: {}.'.format(input_files,output_files))
	cmd =  (
			"{x}: {i}\n" +
			"\t$(ASSEMBLER) -k $(K) \\\n"
			"\t\t-i {ii}\\\n"
			"\t\t-o {oo}\\\n"
			"\t\t-x $@\n\n"
		).format(
			i=' '.join(input_files),
			o=' '.join(output_files),
			ii=' -i '.join(input_files),
			oo=' -o '.join(output_files),
			x=intersection_file,
		)
	print(cmd)


class TreeIndex:

<<<<<<< HEAD
	def __init__(self,tree_newick_fn,index_dir,library_dir,format=1):
=======
	def __init__(self,tree_newick_fn,index_dir,library_dir,format=DEFAULT_FORMAT):
>>>>>>> 09f55ae9e298c875addf2a27dee8e85a33e8da6b
		self.tree_newick_fn=tree_newick_fn
		self.tree=Tree(tree_newick_fn,format=format)
		self.newick_dir=os.path.dirname(tree_newick_fn)
		self.index_dir=index_dir
		self.library_dir=library_dir
		os.makedirs(self.index_dir,exist_ok=True)

	@staticmethod
	def _node_debug(node):
		if hasattr(node,"common_name") and node.common_name!="":
			return "{}_{}".format(node.name,node.common_name)
		elif hasattr(node,"sci_name") and node.sci_name!="":
			return "{}_{}".format(node.name,node.sci_name)
		else:
			return "{}".format(node.name)

	def nonreduced_fasta_fn(self,node):
		return os.path.join(self.index_dir,node.name+".full.fa")
		
	def reduced_fasta_fn(self,node):
		return os.path.join(self.index_dir,node.name+".reduced.fa")

	def process_node(self,node):

		if node.is_leaf():

			if hasattr(node,"fastapath"):
				fastas_fn=node.fastapath.split("@")
				for i in range(len(fastas_fn)):
					fastas_fn[i]=os.path.join(self.library_dir,fastas_fn[i])
				merge_fasta_files(fastas_fn,self.nonreduced_fasta_fn(node),is_leaf=True)

		else:
			#logger.info('BEGIN process non-leaf node "{}"'.format(self._node_debug(node)))
			children=node.get_children()

			# 1) process children
			for child in children:
				self.process_node(child)

			# 2) k-mer propagation & assembly
			input_files=[self.nonreduced_fasta_fn(x) for x in children]
			output_files=[self.reduced_fasta_fn(x) for x in children]
			intersection_file=self.nonreduced_fasta_fn(node)

			# 2a) 1 child
			if len(input_files)==1:
				
				merge_fasta_files(input_files,intersection_file,is_leaf=False)
				#print("ahoj")
				#print(
				#	(
				#		"{new}: {old}\n" +
				#		"\t@cp {old} {new}\n" +
				#		"\n"
				#	).format(new=intersection_file,old=input_files[0])
				#)
				#print(.format(input_files[0],intersection_file))
				#shutil.copyfile(input_files[0],intersection_file)
				#open(output_files[0], 'w').close()

			# 2b) several children
			else:
				assembly(input_files,output_files,intersection_file)
		

	def build_index(self,k,mask_repeats):
		print()
		print("ASSEMBLER=../../bin/assembler")
		print("DUSTMASKER=dustmasker")
		print("K={}".format(k))
		if mask_repeats:
			print("MASK_REPEATS=1")
		print()
		print("ifdef MASK_REPEATS")
		print("   MASKING=$(DUSTMASKER) -infmt fasta -outfmt fasta | sed '/^>/! s/[^AGCT]/N/g'")
		print("else")
		print("   MASKING=tee")
		print("endif")
		print("")
		print("")
		print("all: {}".format(self.nonreduced_fasta_fn(self.tree.get_tree_root())))
		print()
		#logger.info('Going to build index for k={}'.format(k))
		self.process_node(self.tree.get_tree_root())


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Build index.')
	parser.add_argument(
			'-n','--newick-tree',
			type=str,
			metavar='str',
			dest='newick_fn',
			required=True,
			help='Taxonomic tree (Newick).',
		)
	parser.add_argument(
			'-k',
			type=int,
			metavar='int',
			dest='k',
			required=True,
			help='K-mer length k.',
		)
	parser.add_argument(
			'-o','--output-dir',
			type=str,
			metavar='str',
			dest='output_dir_fn',
			required=True,
			help='Output directory (for index).',
		)
	parser.add_argument(
			'-l','--library-dir',
			type=str,
			metavar='str',
			dest='library_dir_fn',
			required=True,
			help='Directory with the library.',
		)
	parser.add_argument(
			'-r','--mask-repeats',
			action='store_true',
			dest='r',
			help='Mask repeats.',
		)

	args = parser.parse_args()

	k=args.k
	assert k>0
	newick_fn=args.newick_fn
	output_dir_fn=args.output_dir_fn
	library_dir_fn=args.library_dir_fn
	r=args.r

	#logger.info("Starting index construction")
	#logger.info("       newick : {}".format(newick_fn))
	#logger.info("   output dir : {}".format(output_dir_fn))
	#logger.info("            k : {}".format(k))

	ti=TreeIndex(
			tree_newick_fn=newick_fn,
			library_dir=library_dir_fn,
			index_dir=output_dir_fn,
		)
	ti.build_index(
			k=k,
			mask_repeats=r,
		)
