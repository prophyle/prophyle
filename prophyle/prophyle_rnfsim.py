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


def parse_args():
    parser = argparse.ArgumentParser(description='Simulator for similarity matrix computation based on RNFtools')

    parser.add_argument(
        'tree_fn',
        type=str,
        metavar='<tree.nhx>',
        help='phylogenetic tree (Newick/NHX)',
    )

    parser.add_argument(
        '-n',
        type=int,
        default=1,
        dest='n_reads',
        help='min number of reads to simulate',
    )

    parser.add_argument(
        '-s',
        choices=simulators,
        default='DwgSim',
        dest='simulator',
        help='simulator to use',
    )

    parser.add_argument(
        '-P',
        action='store_true',
        dest='paired_end',
        help='simulate paired_end reads',
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
