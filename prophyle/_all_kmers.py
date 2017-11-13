#! /usr/bin/env python3
"""Generate all k-mers of a given length.

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT
"""

import argparse
import itertools

parser = argparse.ArgumentParser(description='Generate all possible k-mers of given length.')

parser.add_argument(
    '-k',
    type=int,
    required=True,
    help='k-mer length',
)

parser.add_argument(
    '-f',
    '--format',
    choices=['txt', 'fq', 'fa'],
    default='txt',
    help='output format',
)

args = parser.parse_args()


def print_fq(i, kmer):
    print('@seq' + str(i))
    print(kmer)
    print('+')
    print(args.k * '~')


def print_fa(i, kmer):
    print('>' + str(i))
    print(kmer)


def print_txt(i, kmer):
    print(kmer)


if args.format == "txt":
    pr = print_txt

elif args.format == "fq":
    pr = print_fq

elif args.format == "fa":
    pr = print_fa

i = 1
for kmer in itertools.product('ACGT', repeat=args.k):
    kmer = ''.join(kmer)
    pr(i, kmer)
    i += 1
