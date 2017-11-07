#! /usr/bin/env python3

import argparse
import sys
import os

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro

script_dir = os.path.dirname(os.path.realpath(__file__))
bwa = os.path.join(script_dir, "prophyle_index", "bwa", "bwa")
prophyle_index = os.path.join(script_dir, "prophyle_index", "prophyle_index")


def create_bwa_index(fa):
    # cmd('"{bwa}" index "{fa}"'.format(bwa=bwa,fa=fa))
    pro.run_safe([bwa, 'fa2pac', fa, fa])
    pro.run_safe([bwa, 'pac2bwtgen', fa + ".pac", fa + ".bwt", ">", "/dev/null"])
    pro.run_safe([bwa, 'bwtupdate', fa + ".bwt"])
    pro.run_safe([bwa, 'bwt2sa', fa + ".bwt", fa + ".sa"])


def create_klcp(fa, k):
    pro.run_safe([prophyle_index, 'build', '-k', k, fa, ">", "/dev/null"])


def query(fa, fq, k, u=False, v=False, t=1):
    params = ""
    cmd = [prophyle_index, 'query', "-v" if v else "", "-u" if u else "", '-k', k, '-t', t, fa, fq]
    pro.run_safe(cmd)


def main():
    parser = argparse.ArgumentParser(description='Single-command prophyle_index matching.')
    parser.add_argument(
        '-k',
        type=int,
        metavar='int',
        dest='k',
        required=True,
        help='k-mer length',
    )
    parser.add_argument(
        '-t',
        type=int,
        default=1,
        metavar='int',
        dest='t',
        required=False,
        help='number of threads',
    )
    parser.add_argument(
        '-v',
        action='store_true',
        help='verbose output format',
    )
    parser.add_argument(
        '-u',
        action='store_true',
        help='use rolling window',
    )
    parser.add_argument(
        'in_fasta',
        type=str,
        help='Input FASTA reference.',
    )
    parser.add_argument(
        'in_fq',
        type=str,
        help='Reads to be matched.',
    )

    args = parser.parse_args()

    fa = args.in_fasta
    fq = args.in_fq
    k = args.k
    u = args.u
    v = args.v
    t = args.t

    create_bwa_index(fa)

    if u:
        create_klcp(fa, k)

    query(fa, fq, k, u=u, v=v, t=t)


if __name__ == "__main__":
    main()
