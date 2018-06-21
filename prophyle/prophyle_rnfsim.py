#! /usr/bin/env python3
"""ProPhyle simulator based on RNFtools (https://github.com/karel-brinda/rnftools)

Author:  Simone Pignotti <pignottisimone@gmail.com>

License: MIT

"""

###############################################################################################
###############################################################################################

CONFIG = {
    # read simulator
    'SIMULATOR': 'dwgsim',
    # print diagnostics messages
    'DIAGNOSTICS': False,
}

###############################################################################################
###############################################################################################

import os
import sys
import json
import argparse
import multiprocessing
from ete3 import Tree

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro
import version

simulators = [
    'artillumina',
    'curesim',
    'dwgsim',
    'mason',
    'wgsim',
]

DEFAULT_SIMULATOR = 'dwgsim'
DEFAULT_THREADS = multiprocessing.cpu_count()


def gen_snakefile(fastas, work_dir, simulator, coverage, n_reads, read_len, cluster, paired_end):

    pro.makedirs(work_dir)

    snake_template = os.path.join(os.path.dirname(__file__), 'prophyle_rnfsim_snakefile')
    snake_fn = os.path.join(work_dir, 'Snakefile')
    config_fn = os.path.join(work_dir, 'config.json')
    cluster_fn = os.path.join(work_dir, 'cluster.json')

    if cluster:
        cluster_dict = {
            "__default__": {
                "time": "0-06:00",
                "c": 2,
                "queue": "short",
                "memory": "8G"
            }
        }
        with open(cluster_fn, 'w') as cluster_f:
            json.dump(cluster_dict, cluster_f)

    config_dict = {
        'simulator': simulator,
        'fastas': fastas,
        'reads_in_tuple': 2 if paired_end else 1,
        'read_length': read_len,
        'coverage': coverage,
        'number_of_read_tuples': n_reads,
    }

    with open(config_fn, 'w') as config_f:
        json.dump(config_dict, config_f)

    ln_cmd = ['ln', '-s', '-f', snake_template, snake_fn]
    pro.run_safe(ln_cmd)

def simulate_all_reads(
            tree_fn,
            work_dir='',
            lib_dir='',
            simulator='dwgsim',
            jobs=DEFAULT_THREADS,
            coverage=1,
            n_reads=1000,
            read_len=100,
            cluster=False,
            paired_end=False,
            run=False,
        ):

    tree_dir = os.path.abspath(os.path.dirname(tree_fn))
    if not lib_dir:
        lib_dir = tree_dir
    if not work_dir:
        work_dir = tree_dir

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

    gen_snakefile(fastas, work_dir, simulator, coverage, n_reads, read_len, cluster, paired_end)

    if run:
        if cluster:
            cluster_fn = os.path.join(work_dir, 'cluster.json')
            snk_cmd = [
                'snakemake',
                '-j', jobs,
                '--cluster-config', cluster_fn,
                '--output-wait', '60',
                --cluster,
                "'sbatch -p {{cluster.queue}} -c {{cluster.c}} -t {{cluster.time}} --mem={{cluster.memory}}'"
            ]
        else:
            snk_cmd = [
                'snakemake',
                '-j', jobs,
            ]

        pro.run_safe(snk_cmd)

def parse_args():
    parser = argparse.ArgumentParser(description='Simulator for similarity matrix computation based on RNFtools')

    parser.add_argument(
        'tree_fn',
        type=str,
        metavar='<tree.nhx>',
        help='phylogenetic tree (Newick/NHX)',
    )

    parser.add_argument(
        '-w',
        type=str,
        default='',
        dest='work_dir',
        metavar='STR',
        help='directory where the simulation Snakefile will be created [dir. of tree]',
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
        default=DEFAULT_SIMULATOR,
        dest='simulator',
        metavar='STR',
        help='simulator to use (artillumina, curesim, dwgsim, mason, wgsim) [dwgsim]',
    )

    parser.add_argument(
        '-j',
        type=int,
        default=DEFAULT_THREADS,
        dest='jobs',
        metavar='INT',
        help='number of snakemake jobs (ignored if -R is not set) [auto ({})]'.format(DEFAULT_THREADS),
    )

    parser.add_argument(
        '-P',
        action='store_true',
        dest='paired_end',
        help='simulate paired_end reads [false]',
    )

    parser.add_argument(
        '-C',
        action='store_true',
        dest='cluster',
        help='create configuration files for cluster submission (SLURM only, other environments need manual adjustment) [false]',
    )

    parser.add_argument(
        '-R',
        action='store_true',
        dest='run',
        help='run snakemake after generating the snakefile [false]',
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
            work_dir=args.work_dir,
            lib_dir=args.lib_dir,
            simulator=args.simulator,
            jobs=args.jobs,
            coverage=args.coverage,
            n_reads=args.n_reads,
            read_len=args.read_len,
            cluster=args.cluster,
            paired_end=args.paired_end,
            run=args.run,
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
