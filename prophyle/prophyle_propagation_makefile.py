#! /usr/bin/env python3

"""Create a Makefile for ProPhyle k-mer propagation.

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT

Example:

	prophyle_propagation_makefile


Propagation parameters (in the Makefile, can be changed through CL):
	* NONPROP: no k-mer propagation (sequences for leaves only)
	* REASM: re-assemble sequences in leaves
	* NONDEL: non-deletative propagation, implies REASM
	* MASKREP: mask repeats in leaves


TODO:
	* Check passing parameters and default paths to programs (e.g., prophyle_assembler).
"""


import os
import shutil
import datetime
import sys
import argparse
import textwrap

from ete3 import Tree

DEFAULT_FORMAT = 1


def size_in_mb(file_fn):
	"""Get file size in megabytes.

	Args:
		file_fn (str): File name.
	"""
	return os.path.getsize(file_fn)/(1024**2)


def merge_fasta_files(input_files,output_file,is_leaf):
	"""Print Makefile lines for merging FASTA files and removing empty lines.

	Args:
		input_files (list of str): List of input files.
		output_file (str): Output file.
		is_leaf (str): Is a leaf (i.e., copying must be done).
	"""

	if is_leaf:
		cmd =  (
				"{o}: {i}\n" +
				"\tcat $^ | $(CMD_MASKING) | $(CMD_REASM) > $@\n\n"
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


def assembly(input_files, output_files, intersection_file, count_file="/dev/null"):
	"""Print Makefile lines for running prophyle_assembler.

	Args:
		input_files (list of str): List of input files.
		output_files (list of str): List of output files.
		intersection_file (str): File with intersection.
		count_file (str): File with count statistics.
	"""

	assert(len(input_files)==len(output_files))
	cmd =  textwrap.dedent("""\
			ifdef NONDEL
			CMD_ASM_OUT_{nid} =
			else
			CMD_ASM_OUT_{nid} = -o {oo}
			endif
			
			ifdef NONPROP
			CMD_ASM_{nid} = touch {x} {o}
			else
			CMD_ASM_{nid} = $(PRG_ASM) -S -k $(K) -i {ii} $(CMD_ASM_OUT_{nid}) -x {x} -s {c}
			endif
			
			{x}: {i}
			\t@echo starting propagation for $@
			\t$(CMD_ASM_{nid})
			""".format(
				i=' '.join(input_files),
				o=' '.join(output_files),
				ii=' -i '.join(input_files),
				oo=' -o '.join(output_files),
				x=intersection_file,
				c=count_file,
				nid=intersection_file,
			)
		)
	print(cmd)


class TreeIndex:
	"""Main class for k-mer propagation.
	"""

	def __init__(self,tree_newick_fn,index_dir,library_dir):
		"""Init the class.

		Args:
			tree_newick_fn (str): Tree file name.
			index_dir (str): Directory of the index.
			library_dir (str): Directory with FASTA files.
		"""
		self.tree_newick_fn=tree_newick_fn
		self.tree=Tree(tree_newick_fn,format=DEFAULT_FORMAT)
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
		"""Get name of the full FASTA file.

		Args:
			node: Node of the tree.
		"""
		return os.path.join(self.index_dir,node.name+".full.fa")

	def reduced_fasta_fn(self,node):
		"""Get name of the reduced FASTA file.

		Args:
			node: Node of the tree.
		"""
		return os.path.join(self.index_dir,node.name+".reduced.fa")

	def count_fn(self,node):
		"""Get FASTA name of the file with k-mer counts.

		Args:
			node: Node of the tree.
		"""
		return os.path.join(self.index_dir,node.name+".count.tsv")

	def process_node(self,node):
		"""Recursive function for treating an individual node of the tree.

		Args:
			node: Node of the tree.
		"""

		if node.is_leaf():

			if hasattr(node,"fastapath"):
				fastas_fn=node.fastapath.split("@")
				for i in range(len(fastas_fn)):
					fastas_fn[i]=os.path.join(self.library_dir,fastas_fn[i])
				merge_fasta_files(fastas_fn,self.nonreduced_fasta_fn(node),is_leaf=True)

		else:
			children=node.get_children()

			# 1) process children
			for child in children:
				self.process_node(child)

			# 2) k-mer propagation & assembly
			input_files=[self.nonreduced_fasta_fn(x) for x in children]
			output_files=[self.reduced_fasta_fn(x) for x in children]
			intersection_file=self.nonreduced_fasta_fn(node)
			count_file=self.count_fn(node)


			# 2a) 1 child
			#if len(input_files)==1:
				#
				#merge_fasta_files(input_files,intersection_file,is_leaf=False)
				##print("ahoj")
				##print(
				##	(
				##		"{new}: {old}\n" +
				##		"\t@cp {old} {new}\n" +
				##		"\n"
				##	).format(new=intersection_file,old=input_files[0])
				##)
				##print(.format(input_files[0],intersection_file))
				##shutil.copyfile(input_files[0],intersection_file)
				##open(output_files[0], 'w').close()

			# 2b) several children
			#else:
			assembly(input_files,output_files,intersection_file,count_file)


	def build_index(self,k,mask_repeats):
		"""Print Makefile for the tree.

		Args:
			k (int): K-mer size.
			mask_repeats (bool): Mask repeats using DustMasker.
		"""


		print(textwrap.dedent("""\
				include params.mk\n
				
				PRG_ASM=prophyle_assembler
				PRG_DUST=dustmasker
				
				$(info )
				$(info /------------------------------------------------------------------)
				
				ifdef K
				   $(info | K-mer length:           $(K))
				else
				   $(error | K-mer length is not specified)
				endif
				
				$(info | Assembler:              $(PRG_ASM))
				
				$(info | DustMasker:             $(PRG_DUST))
				
				ifdef MASKREP
				$(info | Masking repeats:        On)
				CMD_MASKING=$(PRG_DUST) -infmt fasta -outfmt fasta | sed '/^>/! s/[^AGCT]/N/g'
				else
				$(info | Masking repeats:        Off)
				CMD_MASKING=tee
				endif
				
				ifdef NONPROP
				$(info | K-mer propagation:      Off)
				else
				$(info | K-mer propagation:      On)
				endif
				
				ifdef NONDEL
				$(info | K-mer propagation mode: Non-deletative)
				REASM=1
				else
				$(info | K-mer propagation mode: Deletative)
				endif
				
				ifdef REASM
				$(info | Re-assembling leaves:   On)
				CMD_REASM=$(PRG_ASM) -S -i - -o -
				else
				$(info | Re-assembling leaves:   Off)
				CMD_REASM=tee
				endif
				$(info \------------------------------------------------------------------)
				$(info )

				all: {root_red}

				{root_red}: {root_nonred}
					ln -s $< $@

				""".format(
							root_nonred=self.nonreduced_fasta_fn(self.tree.get_tree_root()),
							root_red=self.reduced_fasta_fn(self.tree.get_tree_root()),
						)
			))
		#print("all: {}".format(
		#    " ".join(
		#        [self.nonreduced_fasta_fn(x) for x in self.tree.traverse()]
		#        )
		#    )
		#)
		print()
		self.process_node(self.tree.get_tree_root())


def main():
	parser = argparse.ArgumentParser(description='Create Makefile for ProPhyle k-mer propagation.')
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

	ti=TreeIndex(
			tree_newick_fn=newick_fn,
			library_dir=library_dir_fn,
			index_dir=output_dir_fn,
		)
	ti.build_index(
			k=k,
			mask_repeats=r,
		)


if __name__ == "__main__":
	main()
