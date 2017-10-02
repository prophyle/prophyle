#! /usr/bin/env python3

"""ProPhyle assignment algorithm (reference implementation).

Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT
"""

import argparse
import itertools
import os
import sys

import bitarray
import ete3

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro
import version

###############################################################################################
###############################################################################################

# this should be longer than any possible read, because of CRAM (non-tested yet)
FAKE_CONTIG_LENGTH = 42424242

###############################################################################################
###############################################################################################

class Assignment:
    """Class for handling a single read.

    Args:
        output_fo (file): Output file object.
        tree (ete3.Tree): Phylogenetic tree.
        kmer_lca (bool): Simulate LCA on the k-mer level.
        tie_lca (bool): If a tie (multiple winning nodes), compute LCA.
        annotate (bool): If taxonomic info present in the tree, annotate the assignments using SAM tags.
        dont_translate_blocks (bool): ?

    Attributes:
        output_fo (file): Output file object.
        tree (ete3.Tree): Phylogenetic tree.
        k (int): k-mer size.
        kmer_lca (bool): Simulate LCA on the k-mer level.
        tie_lca (bool): If a tie (multiple winning nodes), compute LCA.
        annotate (bool): If taxonomic info present in the tree, annotate the assignments using SAM tags.
        dont_translate_blocks (bool): ?
        krakline_parser (KraklineParser): Parser of kraklines.
        hit_masks_dict (dict): Hit masks
        cov_masks_dict (dict): Cov masks.
        ass_dict (dict): Assignment dictionary.
    """

    def __init__(self, output_fo, tree, kmer_lca=False, tie_lca=False, annotate=False, dont_translate_blocks=False):
        self.output_fo = output_fo
        self.tree = tree
        self.k = tree.k
        self.kmer_lca = kmer_lca
        self.annotate = annotate
        self.tie_lca = tie_lca
        self.dont_translate_blocks = dont_translate_blocks

        self.krakline_parser = KraklineParser ()

        self.hit_masks_dict = {}
        self.cov_masks_dict = {}
        self.ass_dict = {}


    def process_read(self, krakline, form, measure):
        """Process one Kraken-like line.

        Args:
            krakline (str): Kraken-like line.
            form (str): Expected output format (sam/kraken).
            measure (str): Measure (h1/c1/h2/c2).
        """

        self.krakline_parser.parse_krakline(krakline)
        self.blocks_to_masks(self.krakline_parser.kmer_blocks, self.kmer_lca)
        self.compute_assignments(measure)
        self.filter_assignments(measure)
        self.print_assignments(form, measure)


    def blocks_to_masks(self, kmer_blocks, kmer_lca):
        """Extract hit and coverage masks from krakline blocks (without propagation) and store them in self.{hit,cov}masks_dict.

        Args:
            kmer_blocks (list): List of assigned k-mers, i.e., list of (list of node_names, count).
            simulate_kmer_lca (bool): Simulate k-mer LCA on the k-mer level.
        """

        readlen = sum([x[1] for x in kmer_blocks])

        hitmask_len = readlen
        covmask_len = readlen + self.k - 1

        # zero masks
        hitmask_empty = self.bitarray_block(hitmask_len, 0, 0)
        covmask_empty = self.bitarray_block(covmask_len, 0, 0)

        # fast copying (bitarray trick)
        hitmasks_dict = collections.defaultdict(lambda: bitarray.bitarray(hitmask_empty))
        covmasks_dict = collections.defaultdict(lambda: bitarray.bitarray(covmask_empty))

        pos = 0

        for (node_names, count) in kmer_blocks:

            # Kraken output format: 0 and A have special meanings, no blocks
            if node_names != ["0"] and node_names != ["A"]:

                hitmask_1block = bitarray_block(hitmask_len, count, pos)
                covmask_1block = bitarray_block(covmask_len, count+self.k-1, pos)

                if kmer_lca:
                    node_names = [self.ti.lca(*node_names)]

                for node_name in node_names:
                    hitmasks_dict[node_name] |= hitmask_1block
                    covmasks_dict[node_name] |= covmask_1block

            pos += count

        self.hitmasks_dict = hitmasks_dict
        self.covmasks_dict = covmasks_dict


    def compute_assignments (self):
        """Compute assignments and their characteristics from hitmasks and store them in self.ass_dict.
        """

        # TODO: what if no assignment exists
        nodenames=self.hit_masks_dict.keys()
        self.ass_dict={
            nodename: self.evaluate_single_assignment(nodename) for nodename in nodenames
        }


    def evaluate_single_assignment(self, nodename):
        """Evaluate a single assignment.

        Args:
            nodename (str): Name of the node to evaluate.

        Returns:
            assignment (dict): Assignment dictionary.
        """

        # TODO: fix tree (it should be tree index)

        #################################
        # A) Start with the current masks
        #################################
        hitmask = bitarray.bitarray(self.hitmasks_dict[nodename]),
        covmask = bitarray.bitarray(self.covmasks_dict[nodename]),

        node = self.tree.nodename_to_node[nodename]

        #################################
        # B) Update from ancestors
        #################################
        ancestors=self.tree.upnodenames[nodename] & nodenames
        for anc_nodename in ancestors:
            hitmask |= hitmasks[anc_nodename]
            covmask |= covmasks[anc_nodename]

        #################################
        # C) Calculate characteristics
        #################################
        hit = asg['hitmask'].count()
        cov = asg['covmask'].count()

        assignment= {
            # 1. hit count
            'hitmask': hitmask,
            'hitcigar': self.cigar_from_bitmask(hitmask),
            'h1': [hit],
            'hf': [hit / (self.qlen - self.k + 1)],
            'h2': [hit / self.tree.nodename_to_kmercount[nodename]],

            # 2. coverage
            'covmask': covmask,
            'covcigar': self.cigar_from_bitmask(covmask),
            'c1': [cov],
            'cf': [cov / self.readlen],
            'c2': [cov / self.tree.nodename_to_kmercount[nodename]],

            # 3. other values
            'ln' : [self.readlen],
        }

        return assignment


    def filter_assignments(self, measure):
        """Find the best assignments.

        rname=None => unassigned

        Args:
            measure (str): Measure (h1/c1/h2/c2).
        """

        self.max_val = -1
        self.max_nodenames = []

        for nodename in self.ass_dict:
            ass=self.ass_dict[nodename]

            if ass[measure] > self.max_val:
                self.max_val = ass[measure]
                self.max_nodenames = [nodename]

            elif asg[measure] == self.max_val:
                self.max_nodenames.append(nodename)

        for i, rname in enumerate(self.max_nodenames):
            asg = self.asgs[nodename]
            asg['ii'] = i + 1
            asg['is'] = len(self.max_nodenames)


    def make_lca_from_winners(self):
        #
        # todo: test
        #
        # multiple winners => compute LCA and set only those values that are known
        if len(winners) > 1:
            first_winner = self.asgs[winners[0]]
            print("winner rec", first_winner, file=sys.stderr)
            tie_solved = True
            lca = self.tree.lca(winners)
            winners = [lca]
            # fix what if this node exists!
            asg = self.asgs[lca] = {}
            asg['covmask'] = None
            asg['hitmask'] = None
            # asg['covcigar'] = None

            # asg['hitcigar'] = None

            if measure == "h1" or measure == "h2":
                asg['h1'] = first_winner['h1']
                asg['h2'] = first_winner['h2']
                asg['hf'] = first_winner['hf']

                asg['c1'] = None
                asg['c2'] = None
                asg['cf'] = None

            elif measure == "c1" or measure == "c2":
                asg['c1'] = first_winner['c1']
                asg['c2'] = first_winner['c2']
                asg['cf'] = first_winner['cf']

                asg['h1'] = None
                asg['h2'] = None
                asg['hf'] = None


    def print_assignments(self, form, measure, winners):
        """Print the best assignments.

        Args:
            form (str): Expected output format (sam/kraken).
            measure (str): Measure (h1/c1/h2/c2).
        """

        tag_is=len(self.max_rnames)
        for tag_ii, nodename in enumerate(self.max_rnames, 1):
            ass = self.ass_dict[nodename]
            if form == "sam":
                # compute cigar
                if ass['covmask'] is None:
                    ass['covcigar'] = None
                else:
                    ass['covcigar'] = self.cigar_from_bitmask(asg['covmask'])

                if ass['hitmask'] is None:
                    ass['hitcigar'] = None
                else:
                    ass['hitcigar'] = self.cigar_from_bitmask(asg['hitmask'])

                suffix_parts=["ii:i:{}\tis:i:{}".format(tag_ii, tag_is)]
                self.annotate:
                    suffix_parts.append(self.tree.nodename_to_samannot[nodename])
                self.print_sam_line(nodename, "\t".join(suffix_parts))
            elif form == "kraken":
                self.print_kraken_line(nodename)


    @staticmethod
    def cigar_from_bitmask(bitmask):
        """Create a CIGAR from a binary mask.

        Args:
            mask (list): Bitmask.

        Return:
            cigar (str): Cigar string (X/= ops).
        """

        c = []
        mask_bin = map(bool, mask)
        runs = itertools.groupby(mask_bin)
        for run in runs:
            c.append(str(len(list(run[1]))))
            c.append('=' if run[0] else 'X')
        return "".join(c)


    def print_sam_line(self, node_name, suffix):
        """Print a single SAM record.

        Args:
            node_name (str): Node name. None if unassigned.
            suffix (str): Suffix with additional tags.
        """

        tags = []
        qname = self.qname

        if node_name:
            flag = 0
            mapq = "60"
            pos = "1"
            cigar = self.asgs[node_name]['covcigar']
        else:
            flag = 4
            node_name = None
            pos = None
            mapq = None
            cigar = None

        columns = [
            qname,  # QNAME
            str(flag),  # FLAG
            node_name if node_name else "*",  # RNAME
            pos if pos else "0",  # POS
            mapq if mapq else "0",  # MAPQ
            cigar if cigar else "*",  # CIGAR
            "*",  # RNEXT
            "0",  # PNEXT
            "0",  # TLEN
            self.seq if self.seq else "*",  # SEQ
            self.qual if self.qual else "*",  # QUAL
        ]

        if node_name is not None:
            asg = self.asgs[node_name]

            for tag in ['h1', 'h2', 'hf', 'c1', 'c2', 'cf', 'ln']:
                for val in asg[tag]:
                    columns.append("{}:i:{}".format(tag, val))

            columns.append("ii:i:{}".format(asg['ii']))
            columns.append("is:i:{}".format(asg['is']))

            if asg['hitcigar']:
                columns.append("hc:Z:{}".format(asg['hitcigar']))

        print("\t".join(columns), suffix, file=self.output_fo, sep="")


    def print_sam_header(self):
        """Print SAM headers.
        """
        print("@PG", "PN:prophyle", "ID:prophyle", "VN:{}".format(version.VERSION), sep="\t", file=self.output_fo)
        print("@HD", "VN:1.5", "SO:unsorted", sep="\t", file=self.output_fo)
        for node in self.tree.tree.traverse("postorder"):
            self.tree.name_dict[node.name] = node

            try:
                ur = "\tUR:{}".format(node.fastapath)
            except:
                ur = ""

            try:
                sp = "\tSP:{}".format(node.sci_name)
            except:
                sp = ""

            try:
                as_ = "\tAS:{}".format(node.gi)
            except:
                as_ = ""

            if node.name != '':
                print("@SQ\tSN:{rname}\tLN:{rlen}{as_}{ur}{sp}".format(
                    rname=node.name,
                    rlen=FAKE_CONTIG_LENGTH,
                    as_=as_,
                    ur=ur,
                    sp=sp,
                ), file=self.output_fo)


    # TODO: deeply check the entire functino,
    def print_kraken_line(self, node_name):
        """Print a single record in Kraken format.

        Args:
            node_name (str): Node name. None if unassigned.
        """

        if node_name is None:
            stat = "U"
            node_name = "0"
        else:
            stat = "C"

        if self.assignment_lca:
            lca_rnames = []
            for [node_names, count] in self.kmer_blocks:
                assert len(rnames) == 1

                # todo: ???????
                if node_names[0] == "A" or node_names[0] == "0" or self.dont_translate_blocks:
                    taxid = node_names[0]
                else:
                    taxid = int(self.tree.taxid_dict[node_names[0]])
                lca_rnames.extend(count * [taxid])
            c = []
            runs = itertools.groupby(lca_rnames)
            for run in runs:
                c.append("{}:{}".format(
                    str(run[0]),
                    len(list(run[1]))
                ))
            pseudokrakenmers = " ".join(c)

            if stat == "C":
                taxid = str(self.tree.taxid_dict[rname])
            else:
                taxid = "0"

            # todo: check (why is the foramt different ... taxid?)
            columns = [stat, self.qname, taxid, str(self.qlen), pseudokrakenmers]
        # columns=[stat,self.qname,rname,str(self.qlen)," ".join([ "{}:{}".format(",".join(x[0]),x[1]) for x in self.kmer_blocks])]
        else:
            columns = [stat, self.qname, node_name, str(self.qlen), self.krakmers]

        print("\t".join(columns), file=self.output_fo)



###############################################################################################
###############################################################################################


class TreeIndex:
    """Class for an indexed Phylogenetic tree.

    Args:
        tree_newick_fn (str): Filename of the phylogenetic tree.
        k (int): K-mer size.

    Attributes:
        tree (ete3.Tree): Minimal subtree of the original phylogenetic tree.
        k (int): K-mer size.
        nodename_to_node (dict): node name => node.
        nodename_to_taxid (dict): node name => taxid (annotations from the tree).
        nodename_to_samannot (dict): node name => string to append in SAM.
        nodename_to_upnodenames (dict): node name => set of node names of ancestors.
        nodename_to_kmercount (dict): nname => number of k-mers (full set).
    """

    def __init__(self, tree_newick_fn, k):
        tree = pro.load_nhx_tree(tree_newick_fn)
        self.tree = pro.minimal_subtree(tree)

        self.k = k

        self.nodename_to_node = {}
        self.nodename_to_kmercount = {}

        self.nodename_to_taxid = {}
        self.nodename_to_samannot = {}

        self.nodename_to_upnodenames = collections.defaultdict(lambda: set())


        for node in self.tree.traverse("postorder"):
            nodename = node.name
            self.nodename_to_node[nodename] = node
            self.nodename_to_kmercount[nodename] = int(node.kmers_full)

            # annotations
            tags_parts = []
            try:
                tags_parts.extend(["\tgi:Z:", node.gi])
            except AttributeError:
                pass

            try:
                tags_parts.extend(["\tti:Z:", node.taxid])
                self.nodename_to_taxid[nodename] = node.taxid
            except AttributeError:
                pass

            try:
                tags_parts.extend(["\tsn:Z:", node.sci_name])
            except AttributeError:
                pass

            try:
                tags_parts.extend(["\tra:Z:", node.rank])
            except AttributeError:
                pass

            self.nodename_to_samannot[nodename] = "".join(tags_parts)

            # set of upper nodes
            while node.up:
                node = node.up
                self.nodename_to_upnodenames[nodename].add(node.name)


    def lca(self, *node_names):
        """Return LCA for a given list of nodes.

        *node_names (list of str): List of node names.
        """
        assert len(node_names) > 0
        if len(node_names) == 1:
            return node_names[0]
        nodes = list(map(lambda x: self.name_dict[x], node_names))
        lca = nodes[0].get_common_ancestor(nodes)

        if lca.is_root() and len(lca.children) == 1:
            lca = lca.children[0]
        assert lca.name != "" #, [x.name for x in lca.children]

        return lca.name


    @staticmethod
    def bitarray_block(alen, blen, pos):
        """Create a bitarray containing a block of one's.

        Args:
            alen (int): Array length.
            blen (int): Block length (within the array).
            pos (int): Position of the block (0-based).

        Return:
            bitarray (bitarray.bitarray)
        """
        return bitarray.bitarray(pos*[False] + blen*[True] + (alen - pos - blen) * [False])


###############################################################################################
###############################################################################################

class KraklineParser():
    """Class for parsing Kraken-like input into a structure.

    Attributes:
        readname (str): Name of the read.
        raedlen (str): Length of the read.
        seq (str): Sequence of nucleotides. None if unknown.
        qual (str): Sequence of qualities. None if unknown.
        kmer_blocks (list of (list of str, int)): Assigned k-mer blocks, list of (nodenames, count).
    """

    def __init__(self):
        self.readname=None
        self.readlen=None
        self.seq=None
        self.qual=None
        self.kmer_blocks=None


    def parse_krakline(self, krakline):
        """Load a krakline to the current object.

        Args:
            krakline (str): Kraken-like line.
        """

        parts = krakline.strip().split("\t")
        self.readname, _, readlen, krakmers = parts[1:5]
        self.readlen = int(readlen)
        if len(parts) == 7:
            self.seq = parts[5]
            self.qual = parts[6]
        else:
            self.seq = None
            self.qual = None

        # list of (count,list of nodes)
        self.kmer_blocks = []

        for blocks in self.krakmers.split(" "):
            (ids, count) = block.split(":")
            count = int(count)
            nodenames = ids.split(",")
            self.kmer_blocks.append((nodenames, count))


    def check_consistency (self, k):
        """Check consistency of the fields loaded from the krakline.

        Args:
            k (int): k-mer length.
        """
        if self.qlen < k:
            return True

        block_len_sum=sum([x[1] for x in self.kmer_blocks])

        if not self.readlen == block_len_sum + k - 1:
            return False

        if self.seq is not None and len(self.seq)!=self.readlen():
            return False

        if self.qual is not None and len(self.qual)!=self.readlen():
            return False

        return True


###############################################################################################
###############################################################################################


def assign_all_reads(
    tree_fn,
    inp_fo,
    lca,
    form,
    k,
    measure,
    annotate,
    tie,
    dont_translate_blocks,
):
    assert form in ["kraken", "sam"]
    assert k > 1

    tree_index = TreeIndex(
        tree_newick_fn=tree_fn,
        k=k,
    )

    assignment = Assignment(
        output=sys.stdout,
        tree_index=tree_index,
        simulate_lca=lca,
        annotate=annotate,
        tie_lca=tie,
        dont_translate_blocks=dont_translate_blocks,
    )

    if form == "sam":
        assignment.print_sam_header()

    for krakline in inp_fo:
        assignment.process_read(krakline, form=form, measure=measure)


def parse_args():
    parser = argparse.ArgumentParser(description='Implementation of assignment algorithm')

    parser.add_argument('newick_fn',
        type=str,
        metavar='<tree.nhx>',
        help='phylogenetic tree (Newick/NHX)',
    )

    parser.add_argument('k',
        type=int,
        metavar='<k>',
        help='k-mer length',
    )

    parser.add_argument('input_file',
        type=argparse.FileType('r'),
        metavar='<assignments.txt>',
        help='assignments in generalized Kraken format',
    )

    parser.add_argument('-f',
        choices=['kraken', 'sam'],
        default='sam',
        dest='format',
        help='format of output [sam]',
    )

    parser.add_argument('-m',
        choices=['h1', 'c1', 'c2', 'h2'],
        default='h1',
        dest='measure',
        help='measure: h1=hit count, c1=coverage, h2=norm.hit count, c2=norm.coverage [h1]',
    )

    parser.add_argument('-A',
        action='store_true',
        dest='annotate',
        help='annotate assignments',
    )

    parser.add_argument('-L',
        action='store_true',
        dest='tie',
        help='use LCA when tie (multiple hits with the same score)',
    )

    parser.add_argument('-X',
        action='store_true',
        dest='lca',
        help='replace k-mer matches by their LCA',
    )

    parser.add_argument('-D',
        action='store_true',
        dest='donttransl',
        help='do not translate blocks from node to tax IDs',
    )

    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    newick_fn = args.newick_fn
    inp_fo = args.input_file
    lca = args.lca
    form = args.format
    k = args.k
    measure = args.measure
    annotate = args.annotate
    tie = args.tie
    d = args.donttransl

    try:
        assign_all_reads(
            tree_fn=newick_fn,
            inp_fo=inp_fo,
            lca=lca,
            form=form,
            k=k,
            measure=measure,
            annotate=annotate,
            tie=tie,
            dont_translate_blocks=d,
        )

    # Karel: I don't remember why I was considering also IOError here
    # except (BrokenPipeError, IOError):
    except BrokenPipeError:
        # pipe error (e.g., when head is used)
        sys.stderr.close()
        sys.stdout.close()
        exit(0)

    except KeyboardInterrupt:
        pro.message("Error: Keyboard interrupt")
        pro.close_log()
        exit(1)

    finally:
        sys.stdout.flush()
        sys.stderr.flush()


if __name__ == "__main__":
    main()
