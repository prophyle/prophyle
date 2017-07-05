#! /usr/bin/env python3

"""Analyze results of ProPhyle's classification

Authors: Karel Brinda <kbrinda@hsph.harvard.edu>, Simone Pignotti <pignottisimone@gmail.com>

Licence: MIT
"""


import argparse
import os
import sys
import operator

from ete3 import Tree
from collections import Counter

IN_FMTS=['sam','kraken']
STATS=['w','u','wl','ul']
KNOWN_RANKS=['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']

def parse_args():

	desc = """\
Program: prophyle_analyze.py

Analyze results of ProPhyle's classification.
Stats:
	w: weighted assignments
	u: unique assignments, non-weighted
	wl: weighted assignments, propagated to leaves
	ul: unique assignments, propagated to leaves
	"""

	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
									description=desc)

	parser.add_argument('tree',
			metavar='<tree.nhx>',
			type=str,
			help='Newick/NHX tree used for classification'
		)

	input_group = parser.add_mutually_exclusive_group(required=True)

	input_group.add_argument('-i',
			metavar='STR',
			type=str,
			dest='asgs_list',
			nargs='+',
			default=None,
			help="""ProPhyle output files in the format specified with the -f
					option (default SAM, use - for stdin or one file for each
					sample)"""
		)


	input_group.add_argument('-g',
			metavar='STR',
			type=str,
			dest='histograms',
			nargs='+',
			default=None,
			help="""Histograms previously computed using prophyle analyze. Merge
					them and compute OTU table from the result (assignment files
					are not required)"""
		)

	parser.add_argument('-f',
			metavar='STR',
			type=str,
			dest='in_format',
			choices=IN_FMTS,
			default=IN_FMTS[0],
			help='Input format of assignments [{}]'.format(IN_FMTS[0])
		)

	parser.add_argument('-s',
			metavar='STR',
			type=str,
			dest='stats',
			choices=STATS,
			default=STATS[0],
			help='Statistics to use for the computation of histograms [{}]'.format(STATS[0])
		)

	parser.add_argument('-o',
			metavar='STR',
			type=str,
			dest='out_prefix',
			required=True,
			help="""Prefix for output files (the complete file name will be
					<prefix>.tsv for histograms and <prefix>_<otu_suffix>.tsv
					for otu tables)"""
		)

	parser.add_argument('-l',
			type=int,
			metavar='INT',
			dest='max_lines',
			default=-1,
			help='Output only the n nodes with the highest score [all]'
		)

	parser.add_argument('-t',
			type=str,
			metavar='STR',
			dest='otu_suffix',
			default='otu',
			help='Suffix for otu table file [otu]'
		)

	parser.add_argument('-N',
			action='store_true',
			dest='ncbi',
			help="""Use NCBI taxonomic information to calculate abundances at
					every rank for otu tables [default: propagate weighted
					assigments to leaves and do not output internal nodes]"""
		)

	args = parser.parse_args()
	return args


def load_asgs(in_fns, in_format):
	"""Load ProPhyle's assignments in a supported format into a dictionary
	{sample:{read_name:reference_list}}

	Args:
		in_fns (list of str): input filename
		in_format (str): input format
	"""
	asgs={}

	for fn in in_fns:
		if fn=='-':
			assert len(in_fns)==1, "No support for multiple files with stdin"
			base_fn='stdin'
			f=sys.stdin
		else:
			base_fn=os.path.basename(fn)
			if '.' in base_fn:
				base_fn='.'.join(base_fn.split('.')[:-1])
			assert base_fn not in asgs.keys(), "Duplicated input filename"
			f=open(fn,'r')

		asgs[base_fn]={}
		current_asgs=asgs[base_fn]

		for line in f:
			if in_format=='kraken':
				res,read_name,read_ref=line.split('\t')[0:3]
				if res.strip()=='U':
					continue
			else:
				if line.startswith('@'):
					continue
				read_name, _, read_ref=line.split('\t')[0:3]
				if read_ref=='*':
					continue
			try:
				current_asgs[read_name.strip()].append(read_ref.strip())
			except KeyError:
				current_asgs[read_name.strip()]=[read_ref.strip()]

		f.close()

	return asgs


def asgs_to_leaves(tree,asgs):
	"""Propagate all assignments to leaves. Assignments to internal nodes are
	propagated to ALL descendant leaves.

	Args:
		tree (ete3.Tree): tree in Newick/NHX  format used for classification.
		asgs (dict of str:(dict of str:(list of str))): assignments loaded using load_asgs.
	"""
	asgs_to_leaves={}

	for fn, asg_dict in asgs.items():
		asgs_to_leaves[fn]={}
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
			asgs_to_leaves[fn][qname]=list(l)

	return asgs_to_leaves


def compute_histogram(tree, asgs, stats):
	"""Compute histograms with each statistics corresponding to a non-`None` file

	Args:
		tree (ete3.Tree): tree in Newick/NHX format used for classification.
		asgs (dict of str:(dict of str:(list of str))): assignments loaded using load_asgs.
		stats (str): statistics to use.
	"""
	hist={}
	# propagate all assignments to leaves
	if stats.endswith('l'):
		asgs=asgs_to_leaves(tree, asgs)
	for fn, asg_dict in asgs.items():
		hist[fn]=Counter()
		for qname, ref in asg_dict.items():
			# weighted score for multiple assignments
			if stats.startswith('w'):
				no_asg=float(len(ref))
				for taxid in ref:
					hist[fn][taxid]+=1/no_asg
			# unique assignments only
			elif len(ref)==1:
				taxid=ref[0]
				hist[in_fn][taxid]+=1

	return hist


# def compute_node_distances(tree, main_node_name):
#     """Compute distance of the main node from all other nodes in the tree.
#     Returns a dictionary of couples, with the key being the node name,
#     the first element being the real distance of the node, and the second
#     one being the topological distance (e.g. number of nodes between the two).
#
#     Args:
#         tree (ete3.Tree): tree in Newick/NHX format used for classification.
#         main_node_name (str): node from which the distances should be computed.
#     """
#     distances={}
#
#     if main_node_name=="merge_root":
#         main_node=tree
#     else:
#         try:
#             main_node=tree & main_node_name
#         except:
#             print("[prophyle_analyze] Error: node name {} not in the tree".format(main_node_name),
#                     file=sys.stderr)
#             raise
#
#     for node in tree.traverse('postorder'):
#         dist1=main_node.get_distance(node)
#         dist2=main_node.get_distance(node, topology_only=True)
#         distances[node.name]=(dist1,dist2)
#
#     return distances


def print_histogram(histogram, out_file, max_lines):
	"""Print a histogram in tsv format with header, each line containing:
	 - node name;
	 - score of the node in the first sample;
	 - score of the node in the second sample;
	 - ...

	Args:
		histogram (dict of str:Counter): histogram computed using compute_histogram.
		out_file (file): output file.
		max_lines (int): maximum number of nodes to output, ordered by score (all if max_lines < 0).
	"""
	# sum the histograms of each sample to sort assignments by score
	merged_histo=sum(histogram.values(), Counter())
	samples=sorted(histogram.keys())
	# ref_nodes=sorted(set(ref for sample,asg_dict in asgs.items()
	#                 for qname,ref_list in asg_dict.items() for ref in ref_list))

	# header
	print ("\t".join(["#OTU ID"]+samples), file=out_file)
	for i, (node_name,w) in enumerate(sorted(merged_histo.items(), key=operator.itemgetter(1), reverse=True)):
		if max_lines > 0 and not i < max_lines:
			break
		sample_scores=[histogram[sample][node_name] for sample in samples]
		sample_scores=["{:.2f}".format(round(w,2)) if isinstance(w, float) else str(w) for w in sample_scores]
		print("\t".join([node_name]+sample_scores), file=out_file)


def load_histo(in_fns, tree):
	"""Load histogram previously computed using prophyle_analyze.py

	Args:
		in_fns (file): input filenames.
		tree (ete3.Tree): tree in Newick/NHX format used for classification.
	"""
	tree_names=set(node.name for node in tree.traverse('postorder'))
	histo={}
	for fn in in_fns:
		with open(fn, 'r') as f:
			samples=list(map(str.strip, f.readline().split('\t')[1:]))
			for sample in samples:
				assert sample not in histo, "Duplicated sample ID"
				histo[sample]=Counter()
			for line_num, line in enumerate(f):
				fields=list(map(str.strip, line.split('\t')))
				assert len(fields)==(len(samples)+1), "Malformed histogram (check fields at line {})".format(line_num+2)
				node_name=fields[0]
				scores=fields[1:]
				if not node_name in tree_names:
					print("""[prophyle_analyze] Error: node name {} found in the
							histogram {} but not in the tree""".format(n1_name,fn),
							file=sys.stderr)
					raise
				for sample,score in zip(samples,scores):
					histo[sample][node_name]=float(score.strip())
	return histo


# def build_complete_ncbi_tree(tree):
#     """Build the taxonomic tree including internal nodes (for rank dependent evaluation)
#
#     Args:
#         tree (fileobject): file name of the reference tree, without internal nodes
#     """
#
#     ncbi = NCBITaxa()
#     try:
#         taxa=[n.taxid for n in tree.traverse('postorder')]
#     except KeyError:
#         taxa=[n.name for n in tree.traverse('postorder')]
#
#     matched_taxa=False
#     while not matched_taxa:
#         try:
#             tax2rank=ncbi.get_rank(taxa)
#             matched_taxa=True
#         except KeyError as e:
#             # if a taxid is not found, try to build the tree ignoring it
#             taxid_not_found=int(e.args[0])
#             taxa.remove(taxid_not_found)
#             if log:
#                 print('[prophyle_otu_table] ERROR: TaxID ' + str(taxid_not_found) +
#                       ' not found in ETE DB (try updating it)', file=log)
#             continue
#
#     complete_tree=ncbi.get_topology(taxa, intermediate_nodes=True)
#
#     # count of descendants at each rank (used for normalization of assignments propagation)
#     ranked_desc_count=Counter()
#     leaves_count=Counter()
#     for node in complete_tree.traverse('postorder'):
#         taxid=str(node.taxid)
#         rank=str2rank.get(desc.rank, Rank.NO_RANK)
#         leaves_count[taxid]+=len(node.get_leaves())
#         if rank!=Rank.NO_RANK:
#             for anc in node.iter_ancestors():
#                 anc_taxid=str(anc.taxid)
#                 ranked_desc_count[rank][anc_taxid]+=1
#
#     return complete_tree, tax2rank, ranked_desc_count, leaves_count

def ncbi_tree_info(tree):
	# count of descendants at each rank (used for normalization of assignments propagation)
	ranked_desc_count={rank:Counter() for rank in KNOWN_RANKS}
	leaves_count=Counter()
	tax2rank={}
	for node in tree.traverse('postorder'):
		if node.name!='merge_root':
			try:
				taxid=str(node.taxid)
				rank=node.rank
			except AttributeError:
				print("[prophyle_analyze] Warning: missing attributes for node {}".format(node.name),file=sys.stderr)
				continue
			tax2rank[taxid]=rank
			leaves_count[taxid]+=len(node.get_leaves())
			if rank in KNOWN_RANKS:
				for anc in node.iter_ancestors():
					if anc.name!="merge_root":
						try:
							anc_taxid=str(anc.taxid)
							ranked_desc_count[rank][anc_taxid]+=1
						except AttributeError:
							continue

	return tax2rank, ranked_desc_count, leaves_count

def compute_otu_table(histogram, tree, ncbi=False):
	otu_table={}
	if ncbi:
		tax2rank,ranked_desc_count,leaves_count=ncbi_tree_info(tree)
	for sample,histo in histogram.items():
		otu_table[sample]=Counter(histo)
		otu_t=otu_table[sample]
		if ncbi:
			for node in tree.traverse('postorder'):
				if node.name!='merge_root':
					try:
						taxid=str(node.taxid)
						rank=node.rank
					except AttributeError:
						continue
				count=histo[taxid]
				if count!=0:
					for anc in node.iter_ancestors():
						try:
							anc_taxid=str(anc.taxid)
							otu_t[anc_taxid]+=count
						except AttributeError:
							continue
					for desc in node.iter_descendants():
						desc_taxid=str(desc.taxid)
						desc_rank=desc.rank
						# propagate weighted assigment count to descendants
						# (divided by no of nodes at each rank, or leaves count)
						if desc_rank in KNOWN_RANKS:
							otu_t[desc_taxid]+=count/float(ranked_desc_count[desc_rank][taxid])
						elif desc.is_leaf():
							otu_t[desc_taxid]+=count/float(leaves_count[taxid])
					# remove unranked internal nodes from the table (impossible
					# to calculate propagation score accurately)
					if (rank not in KNOWN_RANKS) and (not node.is_leaf()):
						del otu_t[taxid]
		else:
			for node_name,count in histo.items():
				try:
					node=tree & node_name
				except:
					print("[prophyle_analyze] Error: node name {} not in the tree".format(node_name),
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

	return otu_table

# def print_otu_table(otu_t, tree, out_fn, max_lines, fmt="hdf5", ncbi=False):
#     samples=sorted(otu_t.keys())
#     otus=sorted(set(otu for sample, v in otu_t.items() for otu in v))
#     data=[]
#     if ncbi:
#         otus_metadata=[{'taxonomy':node.named_lineage.split('|'),'rank':node.rank}
#                         for taxid in otus for node in tree.search_nodes(taxid=taxid)]
#     else:
# 		#check that is not merge_root
#         features = set()
#         for n in tree.traverse():
#             features |= n.features
#         otus_metadata=[{f:getattr(node,f) for f in features}
#                         for node_name in otus for node in tree.search_nodes(name=node_name)]
#     assert len(otus)==len(otus_metadata), 'Metadata list is too long (duplicate nodes?)'
#     for otu in otus:
#         samples_counts=[]
#         for sample in samples:
#             samples_counts.append(otu_t[sample][otu])
#         data.append(samples_counts)
#     np_data=np.array(data)
#     table = Table(np_data, otus, samples, otus_metadata,
#                 type='OTU table', table_id=out_fn, generated_by='ProPhyle',
#                 create_date=str(dt.now().isoformat()))
#
#     if fmt == "hdf5":
#         open_biom = h5py.File
#     else:
#         open_biom = open
#     with open_biom(out_fn, 'w') as out_f:
#         if fmt == "json":
#             table.to_json(table.generated_by, direct_io=out_f)
#         elif fmt == "tsv":
#             out_f.write(table.to_tsv())
#         else:
#             table.to_hdf5(out_f, table.generated_by)

def print_taxonomy(tree, out_fn):

	taxa2info={}
	for node in tree.traverse('preorder'):
		if node.name!='merge_root':
			try:
				taxid=str(node.taxid).strip()
				rank=node.rank.strip()
				sci_name=node.sci_name.strip()
				if rank in KNOWN_RANKS:
					taxa2info[taxid]=(rank,sci_name)
			except AttributeError:
				continue

	with open(out_fn, 'w') as out_f:
		print('#taxid\t'+'\t'.join(KNOWN_RANKS), file=out_f)
		for node in tree.traverse('preorder'):
			if node.name!='merge_root':
				try:
					lineage=map(str.strip,node.lineage.split('|'))
					rank=node.rank.strip()
				except AttributeError:
					continue
				if rank in KNOWN_RANKS:
					info=['NA']*len(KNOWN_RANKS)
					for t in lineage:
						try:
							r,sn=taxa2info[t]
						except KeyError:
							continue
						try:
							i=KNOWN_RANKS.index(r)
							info[i]=sn
						except ValueError:
							continue
					# t is the last element of lineage, e.g. the taxid of the node
					print(t+'\t'+'\t'.join(info), file=out_f)


def main():
	args=parse_args()
	tree=Tree(args.tree, format=1)

	if args.asgs_list is not None:
		asgs=load_asgs(args.asgs_list,args.in_format)
		histogram=compute_histogram(tree, asgs, args.stats)
		with open(args.out_prefix+'.tsv', "w") as f:
			print_histogram(histogram, f, args.max_lines)

	elif args.histograms is not None:
		histogram=load_histo(args.histograms,tree)

	otu_table=compute_otu_table(histogram,tree,args.ncbi)
	with open("{}_{}.tsv".format(args.out_prefix, args.otu_suffix), 'w') as f:
		print_histogram(otu_table, f, args.max_lines)
	if args.ncbi:
		print_taxonomy(tree,args.out_prefix+"_tax.tsv")

if __name__ == "__main__":
	try:
		main()
	except BrokenPipeError:
		# pipe error (e.g., when head is used)
		sys.stderr.close()
		exit(0)
	except KeyboardInterrupt:
		pro.message("Error: Keyboard interrupt")
		pro.close_log()
		exit(1)