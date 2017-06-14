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
            metavar='<sam_f>'
			type=argparse.FileType('r'),
            help='ProPhyle output in SAM/BAM format'
		)

	parser.add_argument('tree',
    		metavar='<tree.nw>',
    		type=str,
    		help='Newick/NHX tree from the directory of the index used for classification'
		)

	parser.add_argument('-n','--main-node',
			type=str,
            metavar='<node_name>'
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
            metavar='<otu_suffix>'
            dest='otu_suffix',
            default=None
            help='Compute OTU tables from the histograms and write them to fn<i>_otu_suffix.biom'
        )

    parser.add_argument('-h','--from-histogram',
            type=argparse.FileType('r'),
            metavar='<histo>',
            dest='histo',
            default=None,
            help='Compute OTU table from existing histogram (do not compute any histogram)'
        )

	for i in range(1,5):
		parser.add_argument('-{}'.format(i),
				type=str,
				dest='fn{}'.format(i),
				help='Stats {}'.format(i),
				default=None,
				required=False,
			)

	args = parser.parse_args()
	return args


def load_asgs(bam_fn):
	asgs={}
	samfile = pysam.AlignmentFile(bam_fn, "rb")
	for read in samfile.fetch(until_eof=True):
		if read.is_unmapped:
			continue
		try:
			asgs[read.qname].append(read.reference_name)
		except KeyError:
			asgs[read.qname]=[read.reference_name]
	return asgs


def asgs_to_leaves(tree,asgs):
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
	hist=collections.Counter()
	for qname in asgs:
		if len(asgs[qname])==1:
			taxid=asgs[qname][0]
			hist[taxid]+=1
	return hist


def compute_histogram_weighted(asgs):
	hist=collections.Counter()
	for qname in asgs:
		no_asg=len(asgs[qname])
		for taxid in asgs[qname]:
			hist[taxid]+=1.0/no_asg
	return hist


def compute_histograms(tree, asgs):
	"""asgs: rname => list of taxids
	"""

	leaf_asgs=asgs_to_leaves(tree, asgs)
	h1=compute_histogram_unique(asgs)
	h2=compute_histogram_weighted(asgs)
	h3=compute_histogram_unique(leaf_asgs)
	h4=compute_histogram_weighted(leaf_asgs)

	return (h1, h2, h3, h4)


def compute_node_distances(tree, n1_name, n2_name):
	dist1=None
	dist2=None

	if n1_name=="merge_root":
		node1=tree
	else:
		try:
			node1=tree & n1_name
		except:
			print("Error node name: {}".format(n1_name), file=sys.stderr)
			raise

	if n2_name=="merge_root":
		node2=tree
	else:
		try:
			node2=tree & n2_name
		except:
			print("Error node name: {}".format(n2_name), file=sys.stderr)
			raise

	dist1=node1.get_distance(node2)
	dist2=node1.get_distance(node2, topology_only=True)

	return (dist1, dist2)


def print_histogram(tree, main_node_name, histogram, file, lines):
	print ("#node", "hits", "dist", "dist_nodes", sep="\t", file=file)
	for i, (node_name,w) in enumerate(sorted(histogram.items(), key=operator.itemgetter(1), reverse=True)):
		if not i < lines:
			break
		dist1,dist2=compute_node_distances(tree, main_node_name, node_name)
		print(node_name, "{:.2f}".format(round(w,2)) if isinstance(w, float) else w, "{:.2f}".format(round(dist1,2)), dist2, sep="\t", file=file)


def main():
	args=parse_args()

	tree = ete3.Tree(args.tree, format=1)
	asgs=load_asgs(args.sam)

	histograms = compute_histograms(tree, asgs)

	for (histogram, fn) in zip(histograms, (args.fn1, args.fn2, args.fn3, args.fn4) ):
		if fn is not None:
			with open(fn, "w+") as f:
				print_histogram(tree, args.main_node_name, histogram, f, lines=args.lines)


if __name__ == "__main__":
	try:
		main()
	except (BrokenPipeError, IOError):
		sys.stderr.close()
exit(0)
