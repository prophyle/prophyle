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
import enum
import ete3
import numpy as np
from biom.table import Table
from gzip import open as gzip_open
from datetime import datetime as dt

try:
	import h5py
	HAVE_H5PY = True
except ImportError:
	HAVE_H5PY = False

in_fmts=['sam','bam','cram','uncompressed_bam','kraken']
out_fmts=['h5py','json','tsv']

class OrderedEnum(enum.Enum):
	def __ge__(self, other):
		if self.__class__ is other.__class__:
			return self.value >= other.value
		return NotImplemented

	def __gt__(self, other):
		if self.__class__ is other.__class__:
			return self.value > other.value
		return NotImplemented

	def __le__(self, other):
		if self.__class__ is other.__class__:
			return self.value <= other.value
		return NotImplemented

	def __lt__(self, other):
		if self.__class__ is other.__class__:
			return self.value < other.value
		return NotImplemented


class Rank(OrderedEnum):
	NO_RANK = 7
	SPECIES = 6
	GENUS = 5
	FAMILY = 4
	ORDER = 3
	CLASS = 2
	PHYLUM = 1
	KINGDOM = 0


str2rank = {
	'species': Rank.SPECIES,
	'genus': Rank.GENUS,
	'family': Rank.FAMILY,
	'order': Rank.ORDER,
	'class': Rank.CLASS,
	'phylum': Rank.PHYLUM,
	'kingdom': Rank.KINGDOM,
	'no rank': Rank.NO_RANK
}


def parse_args():

	desc = """\
Program: prophyle_analyze.py

Analyze results of ProPhyle's classification.
Stats:
	(1) Unique assignments, non-weighted
	(2) Weighted assignments
	(3) Unique assignments, propagated to leaves
	(4) Weighted assignments, propagated to leaves
	"""

	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
									description=desc)

	parser.add_argument('tree',
			metavar='<tree.nw>',
			type=str,
			help='Newick/NHX tree from the directory of the index used for classification'
		)

	in_group = parser.add_mutually_exclusive_group(required=True)

	in_group.add_argument('-i','--asg-files',
			metavar='<asgs_fn>',
			type=str,
			dest='asgs_list',
			nargs='+',
			default=None,
			help="""ProPhyle output files in the format specified with the -f
					option (use - for stdin, or one file for each sample)"""
		)

	in_group.add_argument('-s','--from-histogram',
			type=argparse.FileType('r'),
			metavar='<histo>',
			dest='histograms',
			nargs='+',
			default=None,
			help="""Compute OTU table from existing histograms (do not compute
					any other) and write it to <histo>_<otu_suffix>.biom
					[<histo>_otu.biom if -o/--otu-suffix is not specified]"""
		)

	parser.add_argument('-f','--in-format',
			metavar='<format>',
			type=str,
			dest='in_format',
			choices=in_fmts,
			default=in_fmts[0],
			help='Input format of assignments ['+in_fmts[0]+']'
		)

	parser.add_argument('-l','--number-of-records',
			type=int,
			metavar='n',
			dest='max_lines',
			default=-1,
			help='Output only the n nodes with the highest score [all]'
		)

	parser.add_argument('-o','--otu-suffix',
			type=str,
			metavar='<otu_suffix>',
			dest='otu_suffix',
			default='otu',
			help="""Compute OTU tables for each specified statistics or input
					histogram and write them to <histo_f[i]>_<otu_suffix>.biom [otu]"""
		)

	parser.add_argument('-u','--out-format',
			metavar='<format>',
			type=str,
			dest='out_format',
			choices=out_fmts,
			default=out_fmts[0],
			help='Output format of OTU tables ['+out_fmts[0]+']'
		)

	parser.add_argument('-n','--ncbi',
			action='store_true',
			dest='ncbi',
			help="""Use NCBI taxonomic information to annotate the OTU tables and
					calculate abundances at every rank [default: propagate
					weighted assigments to leaves and do not output internal
					nodes; metadata from all common features in the tree]"""
		)

	for i in range(1,5):
		parser.add_argument('-{}'.format(i),
				type=str,
				metavar='<histo_fn_{}>'.format(i),
				dest='fn{}'.format(i),
				help="""Output file for histogram computed using stats {}
						[do not compute if not specified]""".format(i),
				default=None,
				required=False,
			)

	args = parser.parse_args()
	return args


def open_asg(in_fn, in_format):
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
	else:
		raise ValueError('Unknown format')
	return in_f


def load_asgs(in_fns, in_format):
	"""Load ProPhyle's assignments in a supported format into a dictionary
	{sample:{read_name:reference_list}}

	Args:
		in_fns (list of str): input filename
		in_format (str): input format
	"""
	asgs={}
	if in_fns==['-']:
		in_files=[('stdin',open_asg(sys.stdin,in_format))]
	else:
		in_files=[]
		for fn in in_fns:
			base_fn=os.path.basename(fn)
			if '.' in base_fn:
				in_files.append(('.'.join(base_fn.split('.')[:-1]),open_asg(fn,in_format)))
			else:
				in_files.append((base_fn,open_asg(fn,in_format)))

	assert len(in_files)==len(set([f[0] for f in in_files])), "Duplicated input filename"
	for in_fn, in_f in in_files:
		asgs[in_fn]={}
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
				asgs[in_fn][read_name.strip()].append(read_ref.strip())
			except KeyError:
				asgs[in_fn][read_name.strip()]=[read_ref.strip()]
		in_f.close()

	return asgs


def asgs_to_leaves(tree,asgs):
	"""Propagate all assignments to leaves. Assignments to internal nodes are
	propagated to ALL descendant leaves.

	Args:
		tree (ete3.Tree): tree in Newick/NHX  format used for classification.
		asgs (dict of str:(dict of str:(list of str))): assignments loaded using load_asgs.
	"""
	asgs2={}

	for in_fn, asg_dict in asgs.items():
		asgs2[in_fn]={}
		for qname, ref in asg_dict.items():
			l=[]
			for tax in ref:
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
					for leaf in n:
						l.append(leaf.name)
			asgs2[in_fn][qname]=list(l)

	return asgs2


def compute_histogram_unique(asgs):
	"""Compute histograms taking into account unique assignments only.

	Args:
		asgs (dict of str:(dict of str:(list of str))): assignments loaded using load_asgs, possibly propagated to the leaves.
	"""
	hist={}
	for in_fn, asg_dict in asgs.items():
		hist[in_fn]=collections.Counter()
		for qname, ref in asg_dict.items():
			if len(ref)==1:
				taxid=ref[0]
				hist[in_fn][taxid]+=1
	return hist


def compute_histogram_weighted(asgs):
	"""Compute histograms assigning weighted scores to multiple assignments.

	Args:
	asgs (dict of str:(dict of str:(list of str))): assignments loaded using load_asgs, possibly propagated to the leaves.
	"""
	hist={}
	for in_fn, asg_dict in asgs.items():
		hist[in_fn]=collections.Counter()
		for qname, ref in asg_dict.items():
			no_asg=float(len(ref))
			for taxid in ref:
				hist[in_fn][taxid]+=1/no_asg
	return hist


def compute_histograms(tree, out_files, asgs):
	"""Compute histograms with each statistics corresponding to a non-`None` file

	Args:
		tree (ete3.Tree): tree in Newick/NHX format used for classification.
		out_files (list of str): list of files (one for each available statistics).
		asgs (dict of str:(dict of str:(list of str))): assignments loaded using load_asgs.
	"""
	histograms=[None]*len(out_files)
	leaf_asgs=asgs_to_leaves(tree, asgs)
	for i, out_f in enumerate(out_files):
		if out_f is not None:
			if i==0:
				histograms[i]=compute_histogram_unique(asgs)
			elif i==1:
				histograms[i]=compute_histogram_weighted(asgs)
			elif i==2:
				histograms[i]=compute_histogram_unique(leaf_asgs)
			elif i==3:
				histograms[i]=compute_histogram_weighted(leaf_asgs)

	return tuple(histograms)


# def compute_node_distances(tree, main_node_name):
# 	"""Compute distance of the main node from all other nodes in the tree.
# 	Returns a dictionary of couples, with the key being the node name,
# 	the first element being the real distance of the node, and the second
# 	one being the topological distance (e.g. number of nodes between the two).
#
# 	Args:
# 		tree (ete3.Tree): tree in Newick/NHX format used for classification.
# 		main_node_name (str): node from which the distances should be computed.
# 	"""
# 	distances={}
#
# 	if main_node_name=="merge_root":
# 		main_node=tree
# 	else:
# 		try:
# 			main_node=tree & main_node_name
# 		except:
# 			print("[prophyle_analyze] Error: node name {} not in the tree".format(main_node_name),
# 					file=sys.stderr)
# 			raise
#
# 	for node in tree.traverse('postorder'):
# 		dist1=main_node.get_distance(node)
# 		dist2=main_node.get_distance(node, topology_only=True)
# 		distances[node.name]=(dist1,dist2)
#
# 	return distances


def print_histogram(histogram, out_file, max_lines):
	"""Print a histogram in tsv format with header, each line containing:
	 - node name;
	 - score of the node in the first sample;
	 - score of the node in the second sample;
	 - ...

	Args:
		histogram (dict of str:collections.Counter): histogram computed using compute_histograms.
		out_file (file): output file.
		max_lines (int): maximum number of nodes to output, ordered by score (all if max_lines < 0).
	"""
	# sum the histograms of each sample to sort assignments by score
	merged_histo=sum(histogram.values(), collections.Counter())
	samples=sorted(histogram.keys())
	# ref_nodes=sorted(set(ref for sample,asg_dict in asgs.items()
	# 				for qname,ref_list in asg_dict.items() for ref in ref_list))

	# header
	print ("\t".join(["node_name"]+samples), file=out_file)
	for i, (node_name,w) in enumerate(sorted(merged_histo.items(), key=operator.itemgetter(1), reverse=True)):
		if max_lines > 0 and not i < max_lines:
			break
		sample_scores=[histogram[sample][node_name] for sample in samples]
		sample_scores=["{:.2f}".format(round(w,2)) if isinstance(w, float) else str(w) for w in sample_scores]
		print("\t".join([node_name]+sample_scores), file=out_file)


def load_histo(in_fn, tree):
	"""Load histogram previously computed using prophyle_analyze.py

	Args:
		in_fn (file): input filename.
		tree (ete3.Tree): tree in Newick/NHX format used for classification.
	"""
	tree_names=set(node.name for node in tree.traverse('postorder'))
	histogram={}
	with open(in_fn, 'r') as in_f:
		samples=in_f.readline().split('\t')[1:]
		for sample in samples:
			histogram[sample]=collections.Counter()
		for line in in_f:
			fields=line.split('\t')
			node_name=fields[0].strip()
			scores=fields[1:]
			if not node_name in tree_names:
				print("""[prophyle_analyze] Error: node name {} found in the
						histogram {} but not in the tree""".format(n1_name,in_fn),
						file=sys.stderr)
				raise
			for sample,score in zip(samples,scores):
				histogram[sample][node_name]=float(score.strip())
	return histogram


def build_complete_ncbi_tree(tree):
	"""Build the taxonomic tree including internal nodes (for rank dependent evaluation)

	Args:
		tree (fileobject): file name of the reference tree, without internal nodes
	"""

	ncbi = NCBITaxa()
	try:
		taxa=[n.taxid for n in tree.traverse('postorder')]
	except KeyError:
		taxa=[n.name for n in tree.traverse('postorder')]

	matched_taxa=False
	while not matched_taxa:
		try:
			tax2rank=ncbi.get_rank(taxa)
			matched_taxa=True
		except KeyError as e:
			# if a taxid is not found, try to build the tree ignoring it
			taxid_not_found=int(e.args[0])
			taxa.remove(taxid_not_found)
			if log:
				print('[prophyle_otu_table] ERROR: TaxID ' + str(taxid_not_found) +
					  ' not found in ETE DB (try updating it)', file=log)
			pass

	complete_tree=ncbi.get_topology(taxa, intermediate_nodes=True)

	# count of descendants at each rank (used for normalization of assignments propagation)
	ranked_desc_count=collections.Counter()
	leaves_count=collections.Counter()
	for node in complete_tree.traverse('postorder'):
		taxid=str(node.taxid)
		rank=str2rank.get(desc.rank, Rank.NO_RANK)
		leaves_count[taxid]+=len(node.get_leaves())
		if rank!=Rank.NO_RANK:
			for anc in node.iter_ancestors():
				anc_taxid=str(anc.taxid)
				ranked_desc_count[rank][anc_taxid]+=1

	return complete_tree, tax2rank, ranked_desc_count, leaves_count


def compute_otu_tables(histograms, tree, ncbi=False):
	if ncbi:
		tree,tax2rank,ranked_desc_count,leaves_count=build_complete_ncbi_tree(tree)
	otu_tables=[{} if histo is not None else None for histo in histograms]
	for otu_table,histogram in zip(otu_tables,histograms):
		if histogram is not None:
			for sample,histo in histogram.items():
				otu_table[sample]=collections.Counter(histo)
				otu_t=otu_table[sample]
				if ncbi:
					for node in tree.traverse('postorder'):
						taxid=str(node.taxid)
						rank=str2rank.get(node.rank, Rank.NO_RANK)
						count=histo[taxid]
						if count!=0:
							for anc in node.iter_ancestors():
								anc_taxid=str(anc.taxid)
								otu_t[anc_taxid]+=count
							for desc in node.iter_descendants():
								desc_taxid=str(desc.taxid)
								desc_rank=str2rank.get(desc.rank, Rank.NO_RANK)
								# propagate weighted assigment count to descendants
								# (divided by no of nodes at each rank, or leaves count)
								if desc_rank!=Rank.NO_RANK:
									otu_t[desc_taxid]+=count/float(ranked_desc_count[desc_rank][taxid])
								elif desc.is_leaf():
									otu_t[desc_taxid]+=count/float(leaves_count[taxid])
							# remove unranked internal nodes from the table (impossible
							# to calculate propagation score accurately)
							if rank==Rank.NO_RANK and (not node.is_leaf()):
								del otu_t[taxid]
				else:
					for node_name,count in histo.items():
						try:
							node=tree & node_name
						except:
							print("[prophyle_analyze] Error: node name {} not in the tree".format(main_node_name),
									file=sys.stderr)
							raise
						# propagate assigments to internal nodes weighted by the number
						# of descendant leaves, then remove them from the otu table
						if not node.is_leaf():
							leaves=node.get_leaves()
							prop_count=count/float(len(leaves))
							for leaf in leaves:
								otu_t[leaf.name]+=prop_count
							del otu_t[node_name]

	return tuple(otu_tables)

def print_otu_table(otu_t, tree, out_fn, max_lines, fmt="hdf5", ncbi=False):
	samples=sorted(otu_t.keys())
	otus=sorted(set(otu for sample, v in otu_t.items() for otu in v))
	data=[]
	if ncbi:
		otus_metadata=[{'taxonomy':node.named_lineage.split('|'),'rank':node.rank}
						for taxid in otus for node in tree.search_nodes(taxid=taxid)]
	else:
		features = set()
		for n in tree.traverse():
			features |= n.features
		otus_metadata=[{f:getattr(node,f) for f in features}
						for node_name in otus for node in tree.search_nodes(name=node_name)]
	assert len(otus)==len(otus_metadata), 'Metadata list is too long (duplicate nodes?)'
	for otu in otus:
		samples_counts=[]
		for sample in samples:
			samples_counts.append(otu_t[sample][otu])
		data.append(samples_counts)
	np_data=np.array(data)
	table = Table(np_data, otus, samples, otus_metadata,
				type='OTU table', table_id=out_fn, generated_by='ProPhyle',
				create_date=str(dt.now().isoformat()))

	if fmt == "hdf5":
		open_biom = h5py.File
	else:
		open_biom = open
	with open_biom(out_fn, 'w') as out_f:
		if fmt == "json":
			table.to_json(table.generated_by, direct_io=out_f)
		elif fmt == "tsv":
			out_f.write(table.to_tsv())
		else:
			table.to_hdf5(out_f, table.generated_by)


# def write_biom(biomT, output_fp, fmt="hdf5", gzip=False):
#     """
#     Write the BIOM table to a file.
#     :type biomT: biom.table.Table
#     :param biomT: A BIOM table containing the per-sample OTU counts and metadata
#                   to be written out to file.
#     :type output_fp str
#     :param output_fp: Path to the BIOM-format file that will be written.
#     :type fmt: str
#     :param fmt: One of: hdf5, json, tsv. The BIOM version the table will be
#                 output (2.x, 1.0, 'classic').
#     """
#     opener = open
#     mode = 'w'
#     if gzip and fmt != "hdf5":
#         if not output_fp.endswith(".gz"):
#             output_fp += ".gz"
#         opener = gzip_open
#         mode = 'wt'
#
#     # HDF5 BIOM files are gzipped by default
#     if fmt == "hdf5":
#         opener = h5py.File
#
#     with opener(output_fp, mode) as biom_f:
#         if fmt == "json":
#             biomT.to_json(biomT.generated_by, direct_io=biom_f)
#         elif fmt == "tsv":
#             biom_f.write(biomT.to_tsv())
#         else:
#             biomT.to_hdf5(biom_f, biomT.generated_by)
#
# return output_fp


def main():
	args=parse_args()
	tree=ete3.Tree(args.tree, format=1)

	if args.asgs_list is not None:
		asgs=load_asgs(args.asgs_list,args.in_format)
		out_files=(args.fn1,args.fn2,args.fn3,args.fn4)
		histograms=compute_histograms(tree, out_files, asgs)
		for (histogram, fn) in zip(histograms, out_files):
			if fn is not None:
				with open(fn+'.tsv', "w+") as f:
					print_histogram(histogram, f, args.max_lines)
	elif args.histo is not None:
		histograms=load_histo(args.histo,tree)

	otu_tables=compute_otu_tables(histograms,tree,args.ncbi)
	for (otu_t, fn) in zip(otu_tables, out_files):
		if otu_t is not None:
			of=fn+'_'+args.otu_suffix+'.biom'
			print_otu_table(otu_t, tree, of, args.max_lines, args.out_format, args.ncbi)

if __name__ == "__main__":
	try:
		main()
	# why catch IOError and exit without error???
	except (BrokenPipeError, IOError):
		sys.stderr.close()
		sys.exit(0)
