#! /usr/bin/env python3
"""Estimate abundances from ProPhyle's assignment using Expectation Maximization

Author:  Simone Pignotti <pignottisimone@gmail.com>

License: MIT

"""

###############################################################################################
###############################################################################################

CONFIG = {
    # print diagnostics messages
    'DIAGNOSTICS': False,
}

###############################################################################################
###############################################################################################

import os
import sys
import pysam
import argparse
import numpy as np

from ete3 import Tree

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro
import version


# likelihood as defined in centrifuge paper
def e_step(ass_mat, ab_vec):
    # matrix of (a_k . C_ik)
    ab_ass_matrix = ass_mat * ab_vec
    # denominator array
    read_sum = ab_ass_matrix.sum(axis=1)
    est_count = (ab_ass_matrix.transpose()/read_sum).sum(axis=1)
    return est_count


def m_step(est_count_vec, len_vec):
    est_ab = est_count_vec/len_vec
    est_ab /= est_ab.sum()
    return est_ab


# EM abundance estimation
def em_abund(ass_mat, len_vec, max_iter=1e5):
    cur_ab = ass_mat.sum(axis=0)
    cur_ab = cur_ab/cur_ab.sum()
    prev_ab = np.zeros(len(cur_ab))
    iter_n = 0

    while iter_n < max_iter and (np.abs(prev_ab - cur_ab)).sum() > 1e-10:
        prev_ab = cur_ab
        est_count = e_step(ass_mat, prev_ab)
        cur_ab = m_step(est_count, len_vec)
        iter_n += 1
        print('{} step done'.format(iter_n), file=sys.stderr)

    return cur_ab

def analyse_assignments(ass_fn, nodes2leaves, vec_pos):

    ref_n = len(vec_pos)
    assignments = np.zeros((1,ref_n))

    ass_f, in_format = pro.open_asg(ass_fn)
    prev_read_name = ""
    cur_ref = []

    if in_format == 'kraken':
        for read in ass_f:
            fields = read.split('\t')
            if fields[0] == 'C':
                read_name = fields[2]
                read_ref = fields[3]
                if read_name == prev_read_name:
                    cur_ref.append(read_ref)
                else:
                    for ref in cur_ref:
                        try:
                            for leaf in nodes2leaves[ref]:
                                assignments[-1,vec_pos[leaf]] = 1
                        except KeyError:
                            print('[prophyle_sim_matrix] Warning: assignments to {} ignored because not in the tree. Are you using the right tree/index?'.format(ref), file=sys.stderr)
                    cur_ref = [read_ref]
                    prev_read_name = read_name
                    assignments = np.append(assignments, np.zeros((1,ref_n)), axis = 0)
    else:
        for read in ass_f.fetch(until_eof=True):
            if not read.is_unmapped:

                read_name = read.qname
                read_ref = read.reference_name
                if read_name == prev_read_name:
                    cur_ref.append(read_ref)
                else:
                    for ref in cur_ref:
                        try:
                            for leaf in nodes2leaves[ref]:
                                assignments[-1,vec_pos[leaf]] = 1
                        except KeyError:
                            print('[prophyle_sim_matrix] Warning: assignments to {} ignored because not in the tree. Are you using the right tree/index?'.format(ref), file=sys.stderr)
                    cur_ref = [read_ref]
                    prev_read_name = read_name
                    assignments = np.append(assignments, np.zeros((1,ref_n)), axis = 0)

    # last assignment
    if len(cur_ref) > 0:
        try:
            for ref in cur_ref:
                for leaf in nodes2leaves[ref]:
                    assignments[-1,vec_pos[leaf]] = 1
        except KeyError:
            print('[prophyle_sim_matrix] Warning: assignments to {} ignored because not in the tree. Are you using the right tree/index?'.format(ref), file=sys.stderr)
    else:
        # remove last zero row
        assignments = assignments[:-1]

    ass_f.close()

    return assignments


def estimate_abundances(tree_fn, asg_fn, len_vec_fn, out_fn):

    tree = Tree(tree_fn, format=1)
    leaves = [leaf.name for leaf in tree]
    nodes2leaves = {node.name: {leaf.name for leaf in node} for node in tree.traverse("postorder")}
    vec_pos = {leaf.name: i for i, leaf in enumerate(tree)}

    count_fn = '.'.join(asg_fn.split('.')[:-1]+['npy'])
    if os.path.isfile(count_fn):
        print("Loading counts from existing npy file", file=sys.stderr)
        map_counts = np.load(count_fn)
    else:
        map_counts = analyse_assignments(asg_fn, nodes2leaves, vec_pos)
        np.save(count_fn, map_counts)

    assert len(leaves) == len(map_counts[0]), "Length of mappings different from #leaves...try to remove <asg.npy> and analyse assignments again using the right tree!"

    len_vec = np.load(len_vec_fn)
    assert len(leaves) == len(len_vec), "Size of genome length vector different from #leaves...have you used the right index/tree?"

    ab_vec = em_abund(map_counts, len_vec)

    with open(out_fn, 'w') as out_f:
        for leaf, ab in zip(leaves, ab_vec):
            print(leaf, ab, sep='\t', file=out_f)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Estimate abundances from ProPhyle's assignment using Expectation Maximization"
    )

    parser.add_argument(
        'tree_fn',
        type=str,
        metavar='<tree.nw>',
        help='taxonomic tree (tree.preliminary.nw)'
    )

    parser.add_argument(
        'asg_fn',
        type=str,
        metavar='<pseudo_aln.bam>',
        help='assignments (output of prophyle classify)'
    )

    parser.add_argument(
        'len_vec_fn',
        type=str,
        metavar='<genome_lengths.npy>',
        help='genome lengths'
    )

    parser.add_argument(
        'out_fn',
        type=str,
        metavar='<output_fn>',
        help='output file'
    )

    parser.add_argument(
        '-c',
        dest='config',
        metavar='STR',
        nargs='*',
        type=str,
        default=[],
        help='configuration (a JSON dictionary)',
    )

    args = parser.parse_args()
    return args


def main():

    args = parse_args()

    global CONFIG
    prophyle_conf_string = pro.load_prophyle_conf(CONFIG, args.config)

    try:
        estimate_abundances(
            tree_fn=args.tree_fn,
            asg_fn=args.asg_fn,
            len_vec_fn=args.len_vec_fn,
            out_fn=args.out_fn,
        )

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
        try:
            sys.stdout.flush()
        except BrokenPipeError:
            pass
        finally:
            try:
                sys.stderr.flush()
            except:
                pass


if __name__ == '__main__':
    main()
