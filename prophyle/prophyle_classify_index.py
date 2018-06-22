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
import json
import argparse
import multiprocessing
from ete3 import Tree

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro
import version

DEFAULT_THREADS = multiprocessing.cpu_count()
DEFAULT_FORMAT = 'bam'

FORMATS = ['sam', 'bam', 'kraken']

def gen_snakefile(work_dir, idx_dir, read_suffix, class_options, targets, cluster):

    pro.makedirs(work_dir)

    snake_template = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'prophyle_classify_snakefile')
    snake_fn = os.path.join(work_dir, 'Snakefile')
    config_fn = os.path.join(work_dir, 'config.json')
    cluster_fn = os.path.join(work_dir, 'cluster.json')

    if cluster:
        cluster_dict = {
            "__default__": {
                "time": "0-06:00",
                "c": 4,
                "queue": "short",
                "memory": "32G"
            }
        }
        with open(cluster_fn, 'w') as cluster_f:
            json.dump(cluster_dict, cluster_f)

    config_dict = {
        'idx': idx_dir,
        'read_suffix': read_suffix,
        'class_options': class_options,
        'targets': targets,
    }

    with open(config_fn, 'w') as config_f:
        json.dump(config_dict, config_f)

    ln_cmd = ['ln', '-s', '-f', snake_template, snake_fn]
    pro.run_safe(ln_cmd, silent=True)


def classify_all_reads(
            idx_dir,
            work_dir='',
            lib_dir='',
            paired_end=False,
            class_options='',
            format=DEFAULT_FORMAT,
            jobs=DEFAULT_THREADS,
            cluster=False,
            run=False,
        ):

    idx_parent_dir = os.path.abspath(os.path.dirname(idx_dir))
    idx_dir = os.path.abspath(idx_dir)
    tree_fn = os.path.join(idx_dir, 'tree.preliminary.nw')
    if not lib_dir:
        lib_dir = idx_parent_dir
    else:
        lib_dir = os.path.abspath(lib_dir)
    if not work_dir:
        work_dir = idx_parent_dir
    else:
        work_dir = os.path.abspath(work_dir)

    tree = Tree(tree_fn, format=1)
    targets = []

    for leaf in tree:
        if hasattr(leaf, 'path'):
            path_list = leaf.path
        elif hasattr(leaf, 'fastapath'):
            path_list = leaf.fastapath
        else:
            print('[prophyle_rnfsim] Warning: leaf {} has no path or fastapath attribute'.format(leaf.name), file=sys.stderr)
            continue

        for p in path_list.split('@'):
            path = '.'.join(p.split('.')[:-1]) if '.' in os.path.basename(p) else p
            targets.append("{}.{}".format(os.path.join(lib_dir, path), format))

    read_suffix = ['1.fq', '2.fq'] if paired_end else ['fq']

    gen_snakefile(work_dir, idx_dir, read_suffix, class_options, targets, cluster)

    if run:
        if cluster:
            cluster_fn = os.path.join(work_dir, 'cluster.json')
            snk_cmd = [
                'cd', work_dir,
                '&&'
                'snakemake',
                '-j', jobs,
                '--cluster-config', cluster_fn,
                '--output-wait', '60',
                --cluster,
                "'sbatch -p {{cluster.queue}} -c {{cluster.c}} -t {{cluster.time}} --mem={{cluster.memory}}'"
            ]
        else:
            snk_cmd = [
                'cd', work_dir,
                '&&'
                'snakemake',
                '-j', jobs,
            ]

        pro.run_safe(snk_cmd)


def parse_args():
    parser = argparse.ArgumentParser(description='Classify reads simulated from the reference genomes to their index, in order to compute a similarity matrix for abundance estimation')

    parser.add_argument(
        'idx_dir',
        type=str,
        metavar='<idx_dir>',
        help='index directory',
    )

    parser.add_argument(
        '-w',
        type=str,
        default='',
        dest='work_dir',
        metavar='STR',
        help='directory where the classification Snakefile will be created (should be the same use for simulation) [parent dir. of idx]',
    )

    parser.add_argument(
        '-g',
        default='',
        dest='lib_dir',
        metavar='STR',
        help='directory with the library sequences [parent dir. of idx]',
    )

    parser.add_argument(
        '-f',
        default=DEFAULT_FORMAT,
        choices=FORMATS,
        dest='format',
        metavar='STR',
        help='classification format (sam, bam, kraken) [bam]'
    )

    parser.add_argument(
        '-o',
        dest='class_options',
        default='',
        metavar='STR',
        help='classification options (-m, -L, -K etc.) [""]'
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
        help='paired end reads [false]',
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
        classify_all_reads(
            idx_dir=args.idx_dir,
            work_dir=args.work_dir,
            lib_dir=args.lib_dir,
            paired_end=args.paired_end,
            class_options=args.class_options,
            format=args.format,
            jobs=args.jobs,
            cluster=args.cluster,
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
