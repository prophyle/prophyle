#! /usr/bin/env python3
"""Normalize a FASTA file (in order to compare different FASTA files).

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT
"""

import argparse
import re
import sys

reg_splitting = re.compile("[^ACGT]")

comp_dict = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
}


def reverse_complement_str(dna):
    reverse_complement = "".join([comp_dict[x] for x in dna[::-1]])
    return reverse_complement


def load_fasta(fasta_fn):
    print("Loading {}".format(fasta_fn), file=sys.stderr)

    sd = {}
    name = None
    seq = []
    with open(fasta_fn) as f:
        for x in f:
            x = x.upper().strip()
            if x != "":
                if x[0] == ">":
                    if name != None:
                        sd[name] = "".join(seq)
                        seq = []
                    name = x[1:]
                else:
                    seq.append(x)
    if name != None:
        sd[name] = "".join(seq)

    print("Loaded {} ({} bp)".format(fasta_fn, sum([len(x) for x in sd.values()])), file=sys.stderr)
    return sd


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='Normalize a FASTA file. This script is useful to check if two FASTA files are equivalent.'
    )

    parser.add_argument(
        '-i',
        '--input',
        help='input fasta file',
        required=True,
    )

    args = parser.parse_args()

    fasta_dict = load_fasta(args.input)

    contigs = list(fasta_dict.values())
    contigs_canonical = []
    for c in contigs:
        cf = c.upper()
        assert set(cf).issubset(set(["A", "C", "G", "T"]))
        cr = reverse_complement_str(c)
        contigs_canonical.append(min(cf, cr))

    contigs_canonical.sort()

    i = 0
    for c in contigs_canonical:
        print(">{}".format(i))
        print(c)
        i += 1
