#! /usr/bin/env python3

"""Analyze results of ProPhyle's classification

Authors: Karel Brinda <kbrinda@hsph.harvard.edu>, Simone Pignotti <pignottisimone@gmail.com>

Licence: MIT
"""


import argparse
import os
import sys
import operator
import pysam

from ete3 import Tree
from ete3.coretype.tree import TreeError
from collections import Counter

IN_FMTS=['sam','bam','cram','uncompressed_bam','kraken','histo']
STATS=['w','u','wl','ul']
KNOWN_RANKS=['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
KRAKEN_RANKS={
    'superkingdom':'K',
    'phylum':'P',
    'class':'C',
    'order':'O',
    'family':'F',
    'genus':'G',
    'species':'S'
}
METAPHLAN_RANKS={
    'unclassified':'u__',
    'superkingdom':'k__',
    'phylum':'p__',
    'class':'c__',
    'order':'o__',
    'family':'f__',
    'genus':'g__',
    'species':'s__'
}

class NotNCBIException(Exception):
    pass

def parse_args():

    desc = """\
Program: prophyle_analyze.py

Analyze results of ProPhyle's classification.
Stats:
w: weighted assignments
u: unique assignments (ignore multiple assignments)
wl: weighted assignments, propagated to leaves
ul: unique assignments, propagated to leaves
"""

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                    description=desc)

    parser.add_argument('index_dir',
            metavar='{index_dir, tree.nw}',
            type=str,
            help='Index directory or phylogenetic tree'
        )

    parser.add_argument('out_prefix',
            metavar='<out_prefix>',
            type=str,
            help="""Prefix for output files (the complete file names will be
                    <out_prefix>_rawhits.tsv for the raw hit counts table and
                    <out_prefix>_otu.tsv for the otu table)"""
        )

    parser.add_argument('input_fns',
            metavar='<input_fn>',
            type=str,
            nargs='+',
            default=None,
            help="""ProPhyle output files whose format is chosen with the -f
                    option. Use '-' for stdin or multiple files with the same
                    format (one per sample)"""
        )

    parser.add_argument('-s',
            metavar=STATS,
            type=str,
            dest='stats',
            choices=STATS,
            default=STATS[0],
            help="""Statistics to use for the computation of histograms:
                    w (default) => weighted assignments;
                    u => unique assignments, non-weighted;
                    wl => weighted assignments, propagated to leaves;
                    ul => unique assignments, propagated to leaves."""
        )

    parser.add_argument('-f',
            metavar=IN_FMTS,
            type=str,
            dest='in_format',
            choices=IN_FMTS,
            default=None,
            help="""Input format of assignments [auto]. If 'histo' is selected the
                    program expects hit count histograms (*_rawhits.tsv)
                    previously computed using prophyle analyze, it merges them
                    and compute OTU table from the result (assignment files are
                    not required)"""
        )

    args = parser.parse_args()
    return args


def open_asg(in_fn, in_format):

    # try to detect kraken and histo formats automatically
    if in_format is None:
        if '.' in in_fn:
            f_ext=in_fn.split('.')[-1]
            if f_ext in IN_FMTS:
                in_format=f_ext
        else:
            with open(in_fn, 'rb') as f:
                f_start=f.read(2)
                if f_start==b'C\t' or f_start==b'U\t':
                    in_format='kraken'
                elif f_start==b'#O':
                    in_format='histo'

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
    # no need to load assignments if input is a histogram -> go to load_histo
    elif in_format=='histo':
        in_f=None
    # let pysam assess the format
    else:
        in_f=pysam.AlignmentFile(in_fn)
    return in_f, in_format


def load_asgs(in_fns, in_format):
    """Load ProPhyle's assignments in a supported format into a dictionary
    {sample:{read_name:reference_list}}

    Args:
        in_fns (list of str): input filename
        in_format (str): input format
    """
    asgs={}
    histograms=[]

    for fn in in_fns:
        if fn=='-':
            assert len(in_fns)==1, "No support for multiple files with stdin"
            assert in_format is not None, "Not able to infer format from stdin"
            base_fn='stdin'
            f=sys.stdin
        else:
            base_fn=os.path.basename(fn)
            if '.' in base_fn:
                base_fn='.'.join(base_fn.split('.')[:-1])
            assert base_fn not in asgs.keys(), "Duplicated input filename"
            f, f_fmt=open_asg(fn,in_format)

        try:
            # if histogram, skip load_asgs and go to load_histo
            if in_format=='histo':
                histograms.append(fn)
                continue
            elif in_format=='kraken':
                read_iterator=(read for read in f)
            # pysam AlignmentFile (sam, bam etc.)
            else:
                read_iterator=(read for read in f.fetch(until_eof=True))

            asgs[base_fn]={}
            current_asgs=asgs[base_fn]
            unclassified=0

            for read in read_iterator:
                if in_format=='kraken':
                    res,read_name,read_ref=read.split('\t')[0:3]
                    if res.strip()=='U':
                        unclassified+=1
                        continue
                else:
                    if read.is_unmapped:
                        unclassified+=1
                        continue
                    read_name=read.qname
                    read_ref=read.reference_name
                try:
                    if read_ref!='merge_root':
                        current_asgs[read_name.strip()].append(read_ref.strip())
                except KeyError:
                    current_asgs[read_name.strip()]=[read_ref.strip()]
        finally:
            if base_fn!='stdin':
                f.close()

    return asgs, histograms, unclassified


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
                if tax!="merge_root":
                    try:
                        n=tree & tax
                    except TreeError:
                        print("[prophyle_analyze] Node {} not found in the tree".format(n1_name), file=sys.stderr)
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
    unique_hist={}
    # propagate all assignments to leaves
    if stats.endswith('l'):
        asgs=asgs_to_leaves(tree, asgs)
    for fn, asg_dict in asgs.items():
        hist[fn]=Counter()
        unique_hist[fn]=Counter()
        for qname, ref in asg_dict.items():
            # weighted score for multiple assignments
            no_asg=float(len(ref))
            for taxid in ref:
                hist[fn][taxid]+=1/no_asg
            # unique assignments only
            if len(ref)==1:
                taxid=ref[0]
                unique_hist[fn][taxid]+=1

    return hist, unique_hist


def print_histogram(histogram, out_f, tax_tree=None):
    """Print a histogram in tsv format with header, each line containing:
     - node name;
     - score of the node in the first sample;
     - score of the node in the second sample;
     - ...

    Args:
        histogram (dict of str:Counter): histogram computed using compute_histogram.
        out_f (file): output file.
    """
    taxa2names={}
    taxa2info={}
    leaves=[]
    if tax_tree is not None:
        leaves=[n.name for n in tax_tree.get_leaves()]
        for node in tax_tree.traverse('preorder'):
            if node.name!='merge_root':
                node_name=node.name
                try:
                    taxid=str(node.taxid)
                except AttributeError:
                    taxid=node_name
                try:
                    rank=node.rank
                except:
                    rank='no rank'
                try:
                    sci_name=node.sci_name
                except:
                    sci_name=node_name
                try:
                    lineage=node.lineage.split('|')
                except:
                    lineage=[taxid]
                taxa2names[taxid]=node_name
                taxa2info[node_name]=(taxid,rank,sci_name,tuple(lineage))

    # sum the histograms of each sample to sort assignments by score
    merged_histo=sum(histogram.values(), Counter())
    samples=sorted(histogram.keys())

    header=["#OTU_ID"]+samples
    if tax_tree is not None:
        header+=KNOWN_RANKS
    print (*header, sep='\t', file=out_f)

    for node_name, w in merged_histo.most_common():
        sample_scores=[histogram[sample][node_name] for sample in samples]
        out_line=[node_name]+["{:.2f}".format(round(c,2)) if isinstance(c, float) else str(c) for c in sample_scores]
        if tax_tree is not None:
            info=['NA']*len(KNOWN_RANKS)
            try:
                taxid,rank,sci_name,lineage=taxa2info[node_name]
            except KeyError:
                continue
            for t in lineage[:-1]:
                try:
                    name=taxa2names[t]
                except KeyError:
                    continue
                try:
                    _,r,sn,_=taxa2info[name]
                except KeyError:
                    continue
                try:
                    i=KNOWN_RANKS.index(r)
                    info[i]=sn
                except ValueError:
                    continue
            if info[-1]=='NA' and node_name in leaves:
                info[-1]=sci_name
            out_line+=info
        print(*out_line, sep='\t', file=out_f)


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
            if len(samples)>len(KNOWN_RANKS) and samples[-len(KNOWN_RANKS):]==KNOWN_RANKS:
                samples=samples[:-len(KNOWN_RANKS)]
            fields_no=len(samples)+1
            for sample in samples:
                assert sample not in histo, "Duplicated sample ID"
                histo[sample]=Counter()
            for line_num, line in enumerate(f):
                fields=list(map(str.strip, line.split('\t')[:fields_no]))
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


def ncbi_tree_info(tree):
    ncbi=True
    # count of descendants at each rank (used for normalization of assignments propagation)
    ranked_desc_count={rank:Counter() for rank in KNOWN_RANKS}
    tax2rank={'1':'no rank'}
    for node in tree.traverse('postorder'):
        if node.name!='merge_root':
            try:
                node_name=str(node.name)
                taxid=str(node.taxid)
                if '-' in node_name:
                    if taxid != node_name.split('-')[1]:
                        print("""[prophyle_analyze] the taxid of node {} is different
                                from its name. Switching to leaf-mode
                                analysis""".format(node_name),
                                file=sys.stderr
                            )
                        ncbi=False
                        return None, None
                elif node_name!=taxid:
                    print("""[prophyle_analyze] the taxid of node {} is different
                            from its name. Switching to leaf-mode
                            analysis""".format(node_name),
                            file=sys.stderr
                        )
                    ncbi=False
                    return None, None
                rank=node.rank
            except AttributeError:
                print("""[prophyle_analyze] taxonomic annotations (taxid or
                        rank) missing for node {}. Switching to leaf-mode
                        analysis""".format(node.name),
                        file=sys.stderr
                    )
                ncbi=False
                return None, None
            tax2rank[node_name]=rank
            if rank in KNOWN_RANKS:
                anc_ranks=[]
                for anc in node.iter_ancestors():
                    if anc.name=="merge_root":
                        ranked_desc_count[rank]['1']+=1
                    else:
                        try:
                            anc_name=str(anc.name)
                            anc_rank=anc.rank
                        except AttributeError:
                            print("""[prophyle_analyze] taxonomic annotations (taxid or
                                    rank) missing for node {}. Switching to leaf-mode
                                    analysis""".format(anc.name),
                                    file=sys.stderr
                                )
                            ncbi=False
                            return None, None
                        ranked_desc_count[rank][anc_name]+=1

    return tax2rank, ranked_desc_count

def compute_otu_table(histogram, tree):
    otu_table={}
    tax2rank,ranked_desc_count=ncbi_tree_info(tree)
    ncbi=False if (tax2rank is None or ranked_desc_count is None) else True
    for sample,histo in histogram.items():
        otu_table[sample]=Counter(histo)
        otu_t=otu_table[sample]
        if ncbi:
            for node in tree.traverse('postorder'):
                if node.name!='merge_root':
                    node_name=node.name
                    rank=node.rank
                    count=histo[node_name]
                    if count!=0:
                        for anc in node.iter_ancestors():
                            if anc.name!='merge_root':
                                otu_t[str(anc.name)]+=count
                        # leaves=node.get_leaves()
                        # for desc in node.iter_descendants():
                        #   desc_taxid=str(desc.name)
                        #   desc_rank=desc.rank
                        #   # propagate weighted assigment count to descendants
                        #   # (divided by no of nodes at each rank, or leaves count)
                        #   if desc_rank in KNOWN_RANKS:
                        #       otu_t[desc_taxid]+=count/float(ranked_desc_count[desc_rank][taxid])
                        #   elif desc.is_leaf():
                        #       otu_t[desc_taxid]+=count/float(len(leaves))
                        # # remove unranked internal nodes from the table (impossible
                        # # to calculate propagation score accurately)
                        # if (rank not in KNOWN_RANKS) and (not node.is_leaf()):
                        #   del otu_t[taxid]
        else:
            print("[prophyle_analyze] Missing some taxonomic tree annotations, switching to leaves-mode")
            for node_name,count in histo.items():
                try:
                    node=tree & node_name
                except TreeError:
                    print("[prophyle_analyze] Error: node {} not in the tree".format(node_name),
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

    return otu_table, ncbi

# def print_otu_table(otu_t, tree, out_fn, max_lines, fmt="hdf5", ncbi=False):
#     samples=sorted(otu_t.keys())
#     otus=sorted(set(otu for sample, v in otu_t.items() for otu in v))
#     data=[]
#     if ncbi:
#         otus_metadata=[{'taxonomy':node.named_lineage.split('|'),'rank':node.rank}
#                         for taxid in otus for node in tree.search_nodes(taxid=taxid)]
#     else:
#       #check that is not merge_root
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

# def print_taxonomy(tree, out_fn):
#
#   taxa2info={}
#   for node in tree.traverse('preorder'):
#       if node.name!='merge_root':
#           try:
#               taxid=str(node.name).strip()
#               rank=node.rank.strip()
#               sci_name=node.sci_name.strip()
#               if rank in KNOWN_RANKS:
#                   taxa2info[taxid]=(rank,sci_name)
#           except AttributeError:
#               continue
#
#   with open(out_fn, 'w') as out_f:
#       print('#taxid\t'+'\t'.join(KNOWN_RANKS), file=out_f)
#       for node in tree.traverse('preorder'):
#           if node.name!='merge_root':
#               try:
#                   lineage=map(str.strip,node.lineage.split('|'))
#                   rank=node.rank.strip()
#               except AttributeError:
#                   continue
#               if rank in KNOWN_RANKS:
#                   info=['NA']*len(KNOWN_RANKS)
#                   for t in lineage:
#                       try:
#                           r,sn=taxa2info[t]
#                       except KeyError:
#                           continue
#                       try:
#                           i=KNOWN_RANKS.index(r)
#                           info[i]=sn
#                       except ValueError:
#                           continue
#                   # t is the last element of lineage, e.g. the taxid of the node
#                   print(t+'\t'+'\t'.join(info), file=out_f)

def print_kraken_report(otu_table, histogram, unclassified, tree, out_f):
    merged_histo=sum(histogram.values(), Counter())
    merged_otu=sum(otu_table.values(), Counter())
    tot_count=float(sum(merged_histo.values())+unclassified)
    if tree.name=='merge_root':
        indentation={n.name:0 for n in tree.children}
    else:
        indentation={tree.name:0}
    unclass_line=["{:.2f}".format(round(unclassified*100/tot_count,2)),
                str(unclassified), str(unclassified), 'U', '0', 'unclassified']
    print(*unclass_line, sep='\t', file=out_f)
    for node in tree.traverse('preorder'):
        if node.up is not None and node.name not in indentation:
            indentation[node.name]=indentation[node.up.name]+1
        node_name=node.name
        if node_name!='merge_root':
            try:
                taxid=str(node.taxid)
            except AttributeError:
                taxid=node_name
            try:
                sci_name=node.sci_name
            except AttributeError:
                sci_name=node_name
            try:
                rank=node.rank
            except AttributeError:
                rank='no rank'
            count=merged_otu[node_name]
            if count > 0:
                hits=merged_histo[node_name]
                # Percentage of reads covered by the clade rooted at this taxon
                out_line=["{:.2f}".format(round(count*100/tot_count,2))]
                # Number of reads covered by the clade rooted at this taxon
                out_line.append("{:.2f}".format(round(count,2)) if isinstance(count, float) else str(count))
                # Number of reads assigned directly to this taxon
                out_line.append("{:.2f}".format(round(hits,2)) if isinstance(hits, float) else str(hits))
                # Rank code
                out_line.append(KRAKEN_RANKS.get(rank,'-'))
                # NCBI taxonomy ID
                out_line.append(taxid)
                # indented scientific name
                out_line.append('  '*indentation[node_name] + sci_name)

                print(*out_line, sep='\t', file=out_f)
    return tot_count

def print_metaphlan_report(otu_table, tot_count, tree, out_f):
    samples = list(otu_table.keys())
    header = ['ID'] + samples
    print(*header, sep='\t', file=out_f)
    for node in tree.traverse('preorder'):
        node_name = node.name
        counts = [otu_table[s][node_name] for s in samples]
        if sum(counts) > 0:
            try:
                taxid = str(node.taxid)
            except AttributeError:
                taxid = node_name
            try:
                lineage = node.lineage.split('|')
            except AttributeError:
                lineage = [taxid]
            # taxonomic string (r__name|r__name...)
            tax_string = ''
            for t in lineage[:-1]:
                try:
                    n = tree.search_nodes(taxid=t)[0]
                except (AttributeError, IndexError) as e:
                    continue
                try:
                    rank = n.rank
                except:
                    rank = 'no rank'
                if rank in KNOWN_RANKS:
                    try:
                        sci_name = n.sci_name
                    except AttributeError:
                        sci_name = n.name
                    prefix = METAPHLAN_RANKS.get(rank,'u__')
                    tax_string += prefix + sci_name + '|'
            try:
                rank = node.rank
            except AttributeError:
                rank = 'no rank'
            try:
                sci_name = node.sci_name
            except AttributeError:
                sci_name = node_name
            prefix = METAPHLAN_RANKS.get(rank,'u__')
            # if lineage is not coherent, just use this node's info
            if taxid != lineage[-1]:
                tax_string = prefix + sci_name
            else:
                tax_string += prefix + sci_name
            # Percentage of reads covered by the clade rooted at this taxon
            out_line=[tax_string]+["{:.5f}".format(round(c*100/tot_count,5)) for c in counts]
            print(*out_line, sep='\t', file=out_f)

def print_centrifuge_report(otu_table, histogram, unique_histogram, tree, out_f):
    header = ['#name', 'taxID', 'taxRank', 'kmerCount', 'numReads', 'numUniqueReads', 'abundance']
    print(*header, sep='\t', file=out_f)
    merged_histo=sum(histogram.values(), Counter())
    merged_unique_histo=sum(unique_histogram.values(), Counter()) if unique_histogram is not None else None
    merged_otu=sum(otu_table.values(), Counter())
    for node in tree.traverse('preorder'):
        node_name=node.name
        count=merged_otu[node_name]
        kmer_count=int(node.kmers_full) if hasattr(node,'kmers_full') else 0
        if count > 0:
            hits=merged_histo[node_name]
            unique_hits=merged_unique_histo[node_name] if merged_unique_histo is not None else '0'
            try:
                rank=node.rank
            except AttributeError:
                rank='unclassified'
            try:
                sci_name=node.sci_name
            except AttributeError:
                sci_name=node_name
            # Scientific name & taxid
            out_line=[sci_name, node_name]
            # Taxonomic rank
            if node.is_leaf():
                out_line.append('leaf')
            else:
                out_line.append(rank)
            # K-mer count (full, before propagation)
            out_line.append(str(kmer_count))
            # Number of reads assigned directly to this taxon (with weighted multiple assignments)
            out_line.append("{:.2f}".format(round(hits,2)) if isinstance(hits, float) else str(hits))
            # Number of reads assigned UNIQUELY to this taxon
            out_line.append(str(int(unique_hits)))
            # Assigned reads/K-mer count ratio
            out_line.append("{:.8f}".format(round(hits/kmer_count,8)) if kmer_count!=0 else 0)
            print(*out_line, sep='\t', file=out_f)

def main():
    args=parse_args()
    if os.path.isdir(args.index_dir):
        tree_path=os.path.join(args.index_dir, 'tree.preliminary.nw')
    else:
        tree_path=args.index_dir
    tree=Tree(tree_path, format=1)

    asgs,histo_fns,unclassified=load_asgs(args.input_fns,args.in_format)
    if len(histo_fns)>0:
        assert len(args.input_fns)==len(histo_fns), "Mixed histogram/assignments input formats not supported"
        histogram=load_histo(histo_fns,tree)
        unique_histogram=None
    else:
        histogram,unique_histogram=compute_histogram(tree, asgs, args.stats)
    otu_table,ncbi=compute_otu_table(histogram,tree)
    with open(args.out_prefix+'.rawhits.tsv', 'w') as f:
        if args.stats.startswith('w'):
            print_histogram(histogram, f, tree if ncbi else None)
        elif unique_histogram is not None:
            print_histogram(unique_histogram, f, tree if ncbi else None)
    with open(args.out_prefix+'.otu.tsv', 'w') as f:
        print_histogram(otu_table, f, tree if ncbi else None)
    if ncbi:
        with open(args.out_prefix+'.kraken.tsv', 'w') as f:
            tot_count=print_kraken_report(otu_table, histogram, unclassified, tree, f)
        with open(args.out_prefix+'.metaphlan.tsv', 'w') as f:
            print_metaphlan_report(otu_table, tot_count, tree, f)
        with open(args.out_prefix+'.centrifuge.tsv', 'w') as f:
            print_centrifuge_report(otu_table, histogram, unique_histogram, tree, f)

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
