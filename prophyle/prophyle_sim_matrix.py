#! /usr/bin/env python3
"""Assign reads simulated from the index to the index itself
(needed for similarity matrix computation)

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
import multiprocessing
import concurrent.futures
from ete3 import Tree

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro
import version

DEFAULT_THREADS = multiprocessing.cpu_count()


def analyse_assignments(leaf_idx, ass_fn_list, nodes2leaves, vec_pos):

    assignments = np.zeros(len(vec_pos))

    for ass_fn in ass_fn_list:
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
                                    assignments[vec_pos[leaf]] += 1
                            except KeyError:
                                print('[prophyle_sim_matrix] Warning: assignments to {} ignored because not in the tree. Are you using the right tree/index?'.format(ref), file=sys.stderr)
                        cur_ref = [read_ref]
                        prev_read_name = read_name
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
                                    assignments[vec_pos[leaf]] += 1
                            except KeyError:
                                print('[prophyle_sim_matrix] Warning: assignments to {} ignored because not in the tree. Are you using the right tree/index?'.format(ref), file=sys.stderr)
                        cur_ref = [read_ref]
                        prev_read_name = read_name

        # last assignment
        if len(cur_ref) > 0:
            try:
                for ref in cur_ref:
                    for leaf in nodes2leaves[ref]:
                        assignments[vec_pos[leaf]] += 1
            except KeyError:
                print('[prophyle_sim_matrix] Warning: assignments to {} ignored because not in the tree. Are you using the right tree/index?'.format(ref), file=sys.stderr)

        ass_f.close()

    # if assignments files don't exist, set the corresponding column to
    # 0 0 .. 1[i] .. 0 0
    if sum(assignments) == 0:
        assignments[leaf_idx] = 1

    # normalize
    assignments /= assignments[leaf_idx]

    # check that the right leaf has the highest number of assignments
    for i, a in enumerate(assignments):
        if a > 1:
            print("[prophyle_sim_matrix] Warning: leaf {} has more assignments than leaf {} for the reads simulated from {}; going to reduce its entry to 1".format(i, leaf_idx, leaf_idx), file=sys.stderr)
            assignments[i] = 1.

    return assignments


def compute_sim_matrix(tree_fn, lib_dir, out_fn, jobs):

    tree = Tree(tree_fn, format=1)
    assigned_fns = []

    nodes2leaves = {node.name: {leaf.name for leaf in node} for node in tree.traverse("postorder")}
    vec_pos = {leaf.name: i for i, leaf in enumerate(tree)}

    for leaf in tree:
        if hasattr(leaf, 'path'):
            path = leaf.path
        elif hasattr(leaf, 'fastapath'):
            path = leaf.fastapath
        else:
            print('[prophyle_sim_matrix] Warning: leaf {} has no path or fastapath attribute'.format(leaf.name), file=sys.stderr)
            continue


        complete_paths = [os.path.join(lib_dir, '.'.join(p.split('.')[:-1])) for p in path.split('@')]
        for i, cp in enumerate(complete_paths):
            if os.path.isfile("{}.kraken".format(cp)):
                complete_paths[i] = "{}.kraken".format(cp)
            elif os.path.isfile("{}.sam".format(cp)):
                complete_paths[i] = "{}.sam".format(cp)
            elif os.path.isfile("{}.bam".format(cp)):
                complete_paths[i] = "{}.bam".format(cp)
            else:
                print('[prophyle_sim_matrix] Warnig: reference genome {} has no assignment file associated (please use the same prefix as the reference, and .kraken, .sam or .bam suffix depending on the classification format)'.format(cp), file=sys.stderr)
                complete_paths[i] = None
        assigned_fns.append(list(complete_paths))

    sim_matrix = np.empty(shape=((len(vec_pos), len(vec_pos))))
    completed = 0

    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
        # Start the load operations and mark each future with its URL
        future2simcolumn = {executor.submit(analyse_assignments, i, ass_fn, nodes2leaves, vec_pos): i for i, ass_fn in enumerate(assigned_fns)}
        for future in concurrent.futures.as_completed(future2simcolumn):
            i = future2simcolumn[future]
            try:
                completed += 1
                column_i = future.result()
            except Exception as e:
                print('[prophyle_sim_matrix] Warning: leaf {} generated an exception: {}. Setting its column to 0 0 .. 1[i] .. 0 0'.format(i, e), file=sys.stderr)
                column_i = np.zeros(len(vec_pos))
                column_i[i] = 1
            print('[prophyle_sim_matrix] {}% - Leaf {} acquired'.format(int(completed*100/len(vec_pos)), i), file=sys.stderr)
            sim_matrix[:,i] = column_i

    np.save(out_fn, sim_matrix)

def parse_args():
    parser = argparse.ArgumentParser(description='Classify reads simulated from the reference genomes to their index, in order to compute a similarity matrix for abundance estimation')

    parser.add_argument(
        'tree_fn',
        type=str,
        metavar='<tree.nhx>',
        help='phylogenetic tree (Newick/NHX)',
    )

    parser.add_argument(
        'lib_dir',
        type=str,
        metavar='<lib_dir>',
        help='library directory containing classification result for each reference genome',
    )

    parser.add_argument(
        'out_fn',
        type=str,
        metavar='<out_fn>',
        help='output filename (matrix is stored in npy format from numpy)',
    )

    parser.add_argument(
        '-j',
        type=int,
        default=DEFAULT_THREADS,
        dest='jobs',
        metavar='INT',
        help='number of threads [auto ({})]'.format(DEFAULT_THREADS),
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
        compute_sim_matrix(
            tree_fn=args.tree_fn,
            lib_dir=args.lib_dir,
            out_fn=args.out_fn,
            jobs=args.jobs,
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


if __name__ == "__main__":
    main()
