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
import numpy
import argparse
import multiprocessing
import concurrent.futures
from ete3 import Tree

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro
import version

DEFAULT_THREADS = multiprocessing.cpu_count()


def analyse_assignments(sam_fn, nodes2leaves, vec_pos):

    assignments = np.array(len(leaves))
    sam_f = pysam.AlignmentFile(sam_fn)
    prev_read_name = ""
    cur_ref = []

    for read in sam_f.fetch(until_eof=True):
        if not read.is_unmapped:
            read_name = read.qname
            read_ref = read.reference_name
            if read_name == prev_read_name:
                cur_ref.append(read_ref)
            else:
                for ref in cur_ref:
                    for leaf in nodes2leaves[ref]:
                        assignments[vec_pos[leaf]] += 1
                cur_ref = [read_ref]
                prev_read_name = read_name

    # last assignment
    if len(cur_ref) > 0:
        for ref in cur_ref:
            for leaf in nodes2leaves[ref]:
                vec_shared[vec_pos[leaf]] += 1

    return assignments


def compute_sim_matrix(tree_fn, lib_dir, out_fn, jobs):

    tree = Tree(tree_fn, format=1)
    assigned_fns = []

    nodes2leaves = {node.name: {leaf.name for leaf in node} for node in tree.traverse("postorder")}
    leaves = [leaf.name for leaf in tree]
    vec_pos = {leaf: i for i, leaf in enumerate(leaves)}

    for leaf in tree:
        if hasattr(leaf, 'path'):
            path = leaf.path
        elif hasattr(leaf, 'fastapath'):
            path = leaf.fastapath
        else:
            print('[prophyle_rnfsim] Warning: leaf {} has no path or fastapath attribute'.format(leaf.name), file=sys.stderr)
            continue

        assert '@' not in path, "[prophyle_sim_matrix] Error: no support for multiple fastas in leaves (triggered by {} in leaf {})".format(path, leaf.name)
        assigned_fns.append(os.path.join(lib_dir, path))

    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
        # Start the load operations and mark each future with its URL
        future_to_simcolumn = {executor.submit(analyse_assignments, ass_fn, nodes2nodes2leaves, vec_pos): ass_fn for ass_fn in assigned_fns}
        for future in concurrent.futures.as_completed(future_to_url):
            url = future_to_url[future]
            try:
                data = future.result()
            except Exception as exc:
                print('%r generated an exception: %s' % (url, exc))
            else:
                print('%r page is %d bytes' % (url, len(data)))

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
            lib_dir=args.idx_dir,
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
