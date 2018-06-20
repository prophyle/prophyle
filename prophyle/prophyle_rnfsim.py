#! /usr/bin/env python3
"""ProPhyle simulator based on RNFtools (https://github.com/karel-brinda/rnftools)

Author:  Simone Pignotti <pignottisimone@gmail.com>

License: MIT

"""

###############################################################################################
###############################################################################################

CONFIG = {
    # read simulator
    'SIMULATOR': 'DwgSim',
    # print diagnostics messages
    'DIAGNOSTICS': False,
}

###############################################################################################
###############################################################################################

import argparse
import os
import sys
from ete3 import Tree

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro
import version

simulators = [
    'ArtIllumina',
    'CuReSim',
    'DwgSim',
    'MasonIllumina',
    'WgSim',
]

def print_snakefile():
    return

def simulate_all_reads(tree_fn,
                        lib_dir='',
                        simulator='DwgSim',
                        coverage=0,
                        n_reads=1000,
                        read_len=100,
                        paired_end=False):
    if not lib_dir:
        lib_dir = os.path.dirname(tree_fn)
    else:
        lib_dir = args.library_dir
    tree = Tree(tree_fn)
    fastas = []
    for leaf in tree:
        if hasattr(leaf, 'path'):
            path_list = leaf.path
        elif hasattr(leaf, 'fastapath'):
            path_list = leaf.fastapath
        else:
            print('[prophyle_rnfsim] Warning: leaf {} has no path or fastapath attribute'.format(leaf.name), file=sys.stderr)
            continue
        for p in path_list.split('@'):
            fastas.append(os.path.join(lib_dir, p))
    print_snakefile(fastas, simulator, coverage, n_reads, read_len, paired_end)

def parse_args():
    parser = argparse.ArgumentParser(description='Simulator for similarity matrix computation based on RNFtools')

    parser.add_argument(
        'tree_fn',
        type=str,
        metavar='<tree.nhx>',
        help='phylogenetic tree (Newick/NHX)',
    )

    parser.add_argument(
        '-g',
        default='',
        dest='lib_dir',
        metavar='STR',
        help='directory with the library sequences [dir. of tree]',
    )

    parser.add_argument(
        '-x',
        type=float,
        default=1.,
        dest='coverage',
        metavar='FLOAT',
        help='min simulation coverage [1.0]',
    )

    parser.add_argument(
        '-n',
        type=int,
        default=1000,
        dest='n_reads',
        metavar='INT',
        help='min number of reads to simulate [1000]',
    )

    parser.add_argument(
        '-l',
        type=int,
        default=100,
        dest='read_len',
        metavar='INT',
        help='reads length [100]',
    )

    parser.add_argument(
        '-s',
        choices=simulators,
        default='DwgSim',
        dest='simulator',
        metavar='STR',
        help='simulator to use (ArtIllumina, CuReSim, DwgSim, MasonIllumina, WgSim) [DwgSim]',
    )

    parser.add_argument(
        '-P',
        action='store_true',
        dest='paired_end',
        help='simulate paired_end reads [false]',
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
        simulate_all_reads(
            tree_fn=args.tree_fn,
            lib_dir=args.lib_dir,
            simulator=args.simulator,
            coverage=args.coverage,
            n_reads=args.n_reads,
            read_len=args.read_len,
            paired_end=args.paired_end,
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
