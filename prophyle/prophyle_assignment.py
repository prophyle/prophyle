#! /usr/bin/env python3

"""ProPhyle assignment algorithm (Python reference implementation).

Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT

TODO:
    * Allow to use c2 and h2 (when k-mer annotations exist in the NHX tree)
"""

import os
import shutil
import datetime
import sys
import argparse
import bitarray
import itertools

import ete3

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro
import version

# this should be longer than any possible read, because of CRAM (non-tested yet)
FAKE_CONTIG_LENGTH = 42424242


def cigar_from_mask(mask):
    """Create a CIGAR from a binary mask.

    Args:
        mask (list): Mask.

    Return:
        cigar (str): Cigar string (X/= ops).
    """

    c = []
    mask_bin=map(bool,mask)
    runs = itertools.groupby(mask_bin)
    for run in runs:
        c.append(str(len(list(run[1]))))
        c.append('=' if run[0] else 'X')
    return "".join(c)


class Read:
    """Class for handling a single read.

    Attributes:
        output (file): Output file object.
        tree (ete3.Tree): Phylogenetic tree.
        simulate_lca (bool): Simulate LCA on the k-mer level.
        tie_lca (bool): If a tie (multiple best nodes), compute LCA.
        annotate (bool): If taxonomic info present in the tree, annotate the assignments using SAM tags.
        dont_translate_blocks (bool): ?
    """


    def __init__(self, output, tree, simulate_lca=False, tie_lca=False, annotate=False, dont_translate_blocks=False):
        self.output = output
        self.tree = tree
        self.k = tree.k
        self.simulate_lca = simulate_lca
        self.annotate = annotate
        self.tie_lca = tie_lca
        self.dont_translate_blocks = dont_translate_blocks


    def process_krakline(self, krakline, form, crit):
        """Process one Kraken-like line.

        Args:
            krakline (str): Kraken-like line.
            form (str): Expected output format (sam/kraken).
            crit (str): Criterion (h1/c1).
        """

        self.load_krakline(krakline)
        self.find_assignments()
        # print(self.asgs)
        self.filter_assignments(crit)
        self.print_assignments(form, crit)


    def load_krakline(self, krakline):
        """Load a krakline to the current object.

        Args:
            krakline (str): Kraken-like line.
        """

        parts = krakline.strip().split("\t")
        self.qname, _, qlen, self.krakmers = parts[1:5]
        self.qlen = int(qlen)
        if len(parts) == 7:
            self.seq = parts[5]
            self.qual = parts[6]
        else:
            self.seq = None
            self.qual = None

        # self.hitmasks=None
        self.asgs = {}

        # list of (count,list of nodes)
        self.kmer_blocks = []

        b_sum = 0
        for b in self.krakmers.split(" "):
            (ids, count) = b.split(":")
            count = int(count)
            b_sum += count
            rnames = ids.split(",")
            if self.simulate_lca:
                rnames = [self.tree.lca(rnames)]
            self.kmer_blocks.append((rnames, count))
        assert self.qlen == b_sum + self.k - 1 or self.qlen < self.k, krakline


    def find_assignments(self):
        """Compute possible assignments of the loaded read.
        """

        # hits before top-down propagation
        hitmasks, covmasks = self.tree.masks_from_kmer_blocks(self.kmer_blocks, simulate_lca=self.simulate_lca)

        rnames = hitmasks.keys()

        # hits after top-down propagation
        for rname in rnames:
            self.asgs[rname] = {
                'hitmask': hitmasks[rname].copy(),
                'covmask': covmasks[rname].copy(),
            }

            node = self.tree.name_dict[rname]
            for p_rname in self.tree.upnodes_dict[rname] & rnames:
                self.asgs[rname]['hitmask'] |= hitmasks[p_rname]
                self.asgs[rname]['covmask'] |= covmasks[p_rname]


    def filter_assignments(self, crit):
        """
        Find the best assignments & compute characteristics (h1, etc.).

        rname=None => unassigned

        Args:
            crit (str): Criterion (h1/c1/...).
        """

        self.max_crit_val = 0
        self.max_crit_rnames = []

        for rname in self.asgs:
            asg = self.asgs[rname]

            """
            1. hit count
            """
            hit = asg['hitmask'].count()
            asg['h1'] = hit
            asg['hf'] = hit/(self.qlen-self.k+1)
            asg['h2'] = hit/self.tree.kmer_count_dict[rname]

            """
            2. coverage
            """
            cov = asg['covmask'].count()
            asg['c1'] = cov
            asg['cf'] = cov/self.qlen
            asg['c2'] = cov/self.tree.kmer_count_dict[rname]

            """
            3. other values
            """
            asg['ln'] = self.qlen

            """
            4. update winners
            """
            if asg[crit] > self.max_crit_val:
                self.max_crit_val = asg[crit]
                self.max_crit_rnames = [rname]
            elif asg[crit] == self.max_crit_val:
                self.max_crit_rnames.append(rname)

        for i,rname in enumerate(self.max_crit_rnames):
            asg = self.asgs[rname]
            asg['ii']=i+1
            asg['is']=len(self.max_crit_rnames)

    def print_assignments(self, form, crit):
        """Print the best assignments.

        Args:
            form (str): Expected output format (sam/kraken).
            crit (str): Criterion (h1/c1).
        """


        winners = self.max_crit_rnames

        # No assignments
        if len(winners)==0:
            if form == "sam":
                self.print_sam_line(None)
            elif form == "kraken":
                self.print_kraken_line(None)
            return

        # Transformation to LCA
        if self.tie_lca:
            #
            # todo: test
            #
            # multiple winners => compute LCA and set only those values that are known
            if len(winners) > 1:
                first_winner=self.asgs[winners[0]]
                tie_solved = True
                lca = self.tree.lca(winners)
                winners = [lca]
                # fix what if this node exists!
                asg = self.asgs[lca] = {}
                asg['covcigar'] = None
                asg['hitcigar'] = None

                if crit == "h1":
                    asg['h1'] = first_winner['h1']
                    asg['h2'] = first_winner['h2']
                    asg['hf'] = first_winner['hf']

                    asg['c1'] = None
                    asg['c2'] = None
                    asg['cf'] = None

                elif crit == "c1":
                    asg['c1'] = first_winner['c1']
                    asg['c2'] = first_winner['c2']
                    asg['cf'] = first_winner['cf']

                    asg['h1'] = None
                    asg['h2'] = None
                    asg['hf'] = None
            # one winner => no
            else:
                pass

        for rname in winners:
            asg = self.asgs[rname]
            if form == "sam":
                # compute cigar
                if self.tie_lca and len(winners) > 1:
                    asg['covcigar'] = None
                    asg['hitcigar'] = None
                else:
                    asg['covcigar'] = cigar_from_mask(asg['covmask'])
                    asg['hitcigar'] = cigar_from_mask(asg['hitmask'])
                self.print_sam_line(rname, self.tree.sam_annotations_dict[rname] if self.annotate else "")
            elif form == "kraken":
                self.print_kraken_line(rname)


    def print_sam_line(self, rname, suffix):
        """Print a single SAM record.

        Args:
            rname (str): Node name. None if unassigned.
            suffix (str): Suffix with additional tags.
        """

        tags = []
        qname = self.qname

        if rname:
            flag = 0
            mapq = "60"
            pos = "1"
            cigar = self.asgs[rname]['covcigar']
        else:
            flag = 4
            rname = None
            pos = None
            mapq = None
            cigar = None

        columns = [
            qname, # QNAME
            str(flag), # FLAG
            rname if rname else "*", # RNAME
            pos if pos else "0", # POS
            mapq if mapq else "0", # MAPQ
            cigar if cigar else "*", # CIGAR
            "*", # RNEXT
            "0", # PNEXT
            "0", # TLEN
            self.seq if self.seq else "*", # SEQ
            self.qual if self.qual else "*", # QUAL
        ]

        if rname is not None:
            asg = self.asgs[rname]

            if asg['h1']:
                columns.append("h1:i:{}".format(asg['h1']))

            if asg['h2']:
                columns.append("h2:f:{}".format(asg['h2']))

            if asg['hf']:
                columns.append("hf:f:{}".format(asg['hf']))

            if asg['c1']:
                columns.append("c1:i:{}".format(asg['c1']))

            if asg['c2']:
                columns.append("c2:f:{}".format(asg['c2']))

            if asg['cf']:
                columns.append("cf:f:{}".format(asg['cf']))


            if asg['ln']:
                columns.append("ln:i:{}".format(asg['ln']))

            if asg['ii']:
                columns.append("ii:i:{}".format(asg['ii']))

            if asg['is']:
                columns.append("is:i:{}".format(asg['is']))

            if asg['hitcigar']:
                columns.append("hc:Z:{}".format(asg['hitcigar']))

        print("\t".join(columns), suffix, file=self.output, sep="")


    def print_sam_header(self):
        """Print SAM headers.
        """

        print("@PG","PN:prophyle","ID:prophyle", "VN:{}".format(version.VERSION), sep="\t", file=self.output)
        print("@HD","VN:1.5","SO:unsorted", sep="\t", file=self.output)
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
                ), file=self.output)


    def print_kraken_line(self, rname):
        """Print a single Kraken-like record.

        Args:
            rname (str): Node name. None is unassigned.
        """

        if rname is None:
            stat = "U"
            rname = "0"
        else:
            stat = "C"

        if self.simulate_lca:
            lca_rnames = []
            for [rnames, count] in self.kmer_blocks:
                assert len(rnames) == 1

                if rnames[0] == "A" or rnames[0] == "0" or self.dont_translate_blocks:
                    taxid = rnames[0]
                else:
                    taxid = int(self.tree.taxid_dict[rnames[0]])
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

            columns = [stat, self.qname, taxid, str(self.qlen), pseudokrakenmers]
        # columns=[stat,self.qname,rname,str(self.qlen)," ".join([ "{}:{}".format(",".join(x[0]),x[1]) for x in self.kmer_blocks])]
        else:
            columns = [stat, self.qname, rname, str(self.qlen), self.krakmers]

        print("\t".join(columns), file=self.output)


class TreeIndex:
    def __init__(self, tree_newick_fn, k):
        self.tree_newick_fn = tree_newick_fn
        tree = pro.load_nhx_tree(tree_newick_fn)
        self.tree = pro.minimal_subtree (tree)

        self.k = k

        self.name_dict = {}

        self.taxid_dict = {}
        self.sam_annotations_dict = {}
        self.upnodes_dict = {}

        self.kmer_count_dict = {}

        for node in self.tree.traverse("postorder"):
            rname = node.name
            self.name_dict[rname] = node
            self.kmer_count_dict[rname] = int(node.kmers_full)

            # annotations
            tags_parts = []
            try:
                tags_parts.extend(["\tgi:Z:", node.gi])
            except AttributeError:
                pass

            try:
                tags_parts.extend(["\tti:Z:", node.taxid])
                self.taxid_dict[rname] = node.taxid
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

            self.sam_annotations_dict[rname] = "".join(tags_parts)

            # set of upper nodes
            self.upnodes_dict[rname] = set()
            while node.up:
                node = node.up
                self.upnodes_dict[rname].add(node.name)

    """
        kmers_assigned_l:
            list of (list_of_nodes, count)

        dict:
            [
                node => hit_vector,
                node => cov_vector
            ]
    """

    def masks_from_kmer_blocks(self, kmers_assigned_l, simulate_lca=False):
        d_h = {}
        d_c = {}

        npos = sum([x[1] for x in kmers_assigned_l])

        h_len = npos
        c_len = npos + self.k - 1

        pos = 0
        for (rname_l, count) in kmers_assigned_l:
            if rname_l != ["0"] and rname_l != ["A"]:
                v_h = bitarray.bitarray(pos * [False] + count * [True] + (npos - pos - count) * [False])
                v_c = bitarray.bitarray(pos * [False] + (count + self.k - 1) * [True] + (npos - pos - count) * [False])

                for rname in rname_l:
                    try:
                        d_h[rname] |= v_h
                        d_c[rname] |= v_c
                    except KeyError:
                        d_h[rname] = v_h.copy()
                        d_c[rname] = v_c.copy()

            pos += count
        return (d_h, d_c)

    def lca(self, noden_l):
        """Return LCA for a given list of nodes.
        """
        assert len(noden_l) > 0, noden_l
        if len(noden_l) == 1:
            return noden_l[0]
        nodes_l = list(map(lambda x: self.name_dict[x], noden_l))
        lca = nodes_l[0].get_common_ancestor(nodes_l)

        if lca.is_root() and len(lca.children) == 1:
            lca = lca.children[0]
        assert lca.name != "", [x.name for x in lca.children]
        return lca.name


def assign(
        tree_fn,
        inp_fo,
        lca,
        form,
        k,
        crit,
        annotate,
        tie,
        dont_translate_blocks,
):
    assert form in ["kraken", "sam"]
    assert k > 1
    ti = TreeIndex(
        tree_newick_fn=tree_fn,
        k=k,
    )

    read = Read(
        output=sys.stdout,
        tree=ti,
        simulate_lca=lca,
        annotate=annotate,
        tie_lca=tie,
        dont_translate_blocks=dont_translate_blocks,
    )
    if form == "sam":
        read.print_sam_header()

    for x in inp_fo:
        read.process_krakline(x, form=form, crit=crit)


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
        choices=['h1', 'c1'],
        default='h1',
        dest='crit',
        help='measure: h1=hitnumber, c1=coverage [h1]',
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
    crit = args.crit
    annotate = args.annotate
    tie = args.tie
    d = args.donttransl

    try:
        assign(
            tree_fn=newick_fn,
            inp_fo=inp_fo,
            lca=lca,
            form=form,
            k=k,
            crit=crit,
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
