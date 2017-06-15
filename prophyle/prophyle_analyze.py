#! /usr/bin/env python3

"""Analyze results of ProPhyle's classification

Authors: Karel Brinda <kbrinda@hsph.harvard.edu>, Simone Pignotti <pignottisimone@gmail.com>

Licence: MIT
"""


import argparse
import os
import sys
import pysam
import biom
import collections
import operator
import ete3


formats=['sam','bam','cram','uncompressed_bam','kraken']

def parse_args():

	desc = """\
		Program: prophyle_analyze.py

		Analyze results of ProPhyle's classification
		Stats:
			(1) Unique assignments, non-weighted
			(2) Weighted assignments
			(3) Unique assignments, propagated to leaves
			(4) Weighted assignments, propagated to leaves
	"""

	parser = argparse.ArgumentParser(description=desc)

	parser.add_argument('sam',
			metavar='<sam_f>',
			type=argparse.FileType('r'),
			help='ProPhyle output in the format specified with the -f option'
		)

	parser.add_argument('tree',
			metavar='<tree.nw>',
			type=str,
			help='Newick/NHX tree from the directory of the index used for classification'
		)

	parser.add_argument('-f','--format',
			type=str,
			metavar='<format>',
			choices=formats,
			default=formats[0]
			help='Input format ['+formats[0]+']'
		)

	parser.add_argument('-n','--main-node',
			type=str,
			metavar='<node_name>',
			dest='main_node_name',
			default='merge_root',
			help='Node from which the distances should be computed [root]'
		)

	parser.add_argument('-l','--number-of-records',
			type=int,
			metavar='n',
			dest='lines',
			default=-1,
			help='Output only the n nodes with the highest score [all]'
		)

	parser.add_argument('-t','--otu',
			type=str,
			metavar='<otu_suffix>',
			dest='otu_suffix',
			default=None,
			help="""Compute OTU tables for each specified statistics and write
					them to <histo_f[i]>_<otu_suffix>.biom"""
		)

	parser.add_argument('-s','--from-histogram',
			type=argparse.FileType('r'),
			metavar='<histo>',
			dest='histo',
			default=None,
			help="""Compute OTU table from existing histogram (do not compute
					any other histogram) and write it to <histo>_<otu_suffix>.biom
					[<histo>_otu.biom]"""
		)

	for i in range(1,5):
		parser.add_argument('-{}'.format(i),
				type=str,
				metavar='<histo_f[{}]>'.format(i),
				dest='fn{}'.format(i),
				help='Stats {}'.format(i),
				default=None,
				required=False,
			)

	args = parser.parse_args()
	return args


def load_asgs(in_fn, in_format):
	"""Load ProPhyle's assignments in a supported format into a dictionary {read_name : reference_list}

	Args:
		in_fn (str): input filename
		in_format (str): input format
	"""
	asgs={}
	if in_format=='sam':
		in_f=pysam.AlignmentFile(in_fn, "r")
	elif in_format=='bam':
		in_f=pysam.AlignmentFile(in_fn, "rb")
	elif in_format=='cram':
		in_f=pysam.AlignmentFile(in_fn, "rc")
	elif in_format=='uncompressed_bam':
		in_f=pysam.AlignmentFile(in_fn, "ru")
	elif in_format=='kraken':
		in_f=open(in_fn,'r')

	if in_format=='kraken':
		read_iterator=(read for read in in_f)
	else:
		read_iterator=(read for read in in_f.fetch(until_eof=True))

	for read in read_iterator:
		if in_format=='kraken':
			res,read_name,read_ref=read.split('\t')[0:3]
			if res.strip()=='U':
				continue
		else:
			if read.is_unmapped:
				continue
			read_name=read.qname
			read_ref=read.reference_name
		try:
			asgs[read_name].append(read_ref)
		except KeyError:
			asgs[read_name]=[read_ref]

	return asgs


def asgs_to_leaves(tree,asgs):
	"""Propagate all assignments to leaves. Assignments to internal nodes are
	propagated to ALL descendant leaves.

	Args:
		tree (ete3.Tree): tree in Newick/NHX  format used for classification.
		asgs (dict of tuple): Assignments loaded using load_asgs.
	"""
	asgs2={}

	for qname in asgs:
		l=[]
		for tax in asgs[qname]:
			if tax=="merge_root":
				n=tree
			else:
				try:
					n=tree & tax
				except:
					print("Error node name: {}".format(n1_name), file=sys.stderr)
					raise
			if n.is_leaf():
				l.append(n.name)
			else:
				for desc in n.iter_descendants("postorder"):
					if desc.is_leaf():
						l.append(desc.name)
			asgs2[qname]=l

	return asgs2


def compute_histogram_unique(asgs):
	"""Compute histograms taking into account unique assignments only.

	Args:
		asgs (dict of tuple): assignments loaded using load_asgs, possibly propagated to the leaves.
	"""
	hist=collections.Counter()
	for qname in asgs:
		if len(asgs[qname])==1:
			taxid=asgs[qname][0]
			hist[taxid]+=1
	return hist


def compute_histogram_weighted(asgs):
	"""Compute histograms assigning weighted scores to multiple assignments.

	Args:
		asgs (dict of tuple): assignments loaded using load_asgs, possibly propagated to the leaves.
	"""
	hist=collections.Counter()
	for qname in asgs:
		no_asg=len(asgs[qname])
		for taxid in asgs[qname]:
			hist[taxid]+=1.0/no_asg
	return hist


def compute_histograms(tree, out_files, asgs):
	"""Compute histograms with each statistics corresponding to a non-`None` file

	Args:
		tree (ete3.Tree): tree in Newick/NHX format used for classification.
		out_files (list of str): list of files (one for each available statistics)
	"""
	histograms=[None]*len(out_files)
	leaf_asgs=asgs_to_leaves(tree, asgs)
	for i, out_f in enumerate(out_files):
		if out_f is not None:
			if i==0:
				histograms[i]=compute_histogram_unique(asgs)
			elif i==1:
				histograms[i]=compute_histogram_unique(asgs)
			elif i==2:
				histograms[i]=compute_histogram_unique(asgs)
			elif i==3:
				histograms[i]=compute_histogram_unique(asgs)

	return tuple(histograms)


def compute_node_distances(tree, main_node_name):
	"""Compute distance of the main node from all other nodes in the tree.
	Returns a dictionary of couples, with the key being the node name,
	the first element being the real distance of the node, and the second
	one being the topological distance (e.g. number of nodes between the two).

	Args:
		tree (ete3.Tree): tree in Newick/NHX format used for classification.
		main_node_name (str): node from which the distances should be computed.
	"""
	distances={}

	if main_node_name=="merge_root":
		main_node=tree
	else:
		try:
			main_node=tree & main_node_name
		except:
			print("[prophyle_analyze] Error: node name {} not in the tree".format(main_node_name),
					file=sys.stderr)
			raise

	for node in tree.iterate('postorder'):
		dist1=main_node.get_distance(node)
		dist2=main_node.get_distance(node, topology_only=True)
		distances[node.name]=(dist1,dist2)


def print_histogram(histogram, distances, out_file, lines):
	"""Print a histogram in tsv format with header, each line containing:
	 - node name;
	 - score of the node;
	 - distance of the node from the main node;
	 - topological distance of the node from the main node.

	Args:
		histogram (dict of float): histogram computed using `compute_histograms`.
		distances (dict of tuple): distance of each node computed using `compute_node_distances`.
		out_file (file): output file.
		lines (int): maximum number of nodes to output, ordered by score (all if lines < 0).
	"""
	# maybe score is better than hits, and something like unit_dist is better than dist_nodes
	print ("node_name", "hits", "dist", "dist_nodes", sep="\t", file=out_file)
	# reverse = true???
	for i, (node_name,w) in enumerate(sorted(histogram.items(), key=operator.itemgetter(1), reverse=True)):
		if i > 0 and not i < lines:
			break
		try:
			dist1,dist2=distances[node_name]
			print(node_name, "{:.2f}".format(round(w,2)) if isinstance(w, float) else w,
					"{:.2f}".format(round(dist1,2)), dist2, sep="\t", file=out_file)
		except KeyError:
			print("[prophyle_analyze] Error: node name {} found in assignments but not in the tree".format(n1_name),
					file=sys.stderr)
			raise


def load_histo(in_fn, distances):
	"""Load histogram previously computed using prophyle_analyze.py

	Args:
		in_fn (file): input filename.
		distances (dict of tuple): distance of each node computed using `compute_node_distances`.
	"""
	histogram=collections.Counter()
	with open(in_fn, 'r') as in_f:
		for line in in_f:
			node_name,w=line.split('\t')[0:2]
			if not node_name in distances:
				print("[prophyle_analyze] Error: node name {} found in the histogram {} but not in the tree".format(n1_name,in_fn),
						file=sys.stderr)
				raise
			histogram[node_name]=w
	return histogram


def compute_otu_tables():
	return

def main():
	args=parse_args()
	tree=ete3.Tree(args.tree, format=1)
	distances=compute_node_distances(tree,args.main_node_name)

	if args.sam:
		asgs=load_asgs(args.sam)
		histograms=compute_histograms(tree, asgs)
		for (histogram, fn) in zip(histograms, (args.fn1, args.fn2, args.fn3, args.fn4) ):
			if fn is not None:
				with open(fn, "w+") as f:
					print_histogram(tree, args.main_node_name, histogram, f, lines=args.lines)
	elif args.histo:
		histograms=load_histo(args.histo,distances)
	otu_tables=compute_otu_tables()

if __name__ == "__main__":
	try:
		main()
	# why catch IOError and exit without error???
	except (BrokenPipeError, IOError):
		sys.stderr.close()
		sys.exit(0)
