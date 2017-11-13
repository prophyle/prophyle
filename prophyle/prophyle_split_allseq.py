#! /usr/bin/env python3
"""Split fasta files containing sequences of different organisms (with different TaxID)
    Used for HMP dataset, since all sequences are distributed in a single fasta

    Author: Simone Pignotti <pignottisimone@gmail.com>

Example:
    $ prophyle_split_allseq.py <output_dir> [-i <input.fa>]

TODO:
    * support infasta_offset and base_len (already included in the trees) for index construction
"""

import os
import sys
import argparse


def split_fs(output_dir_fn):
    i = 0
    while True:
        i += 1
        fn = os.path.join(output_dir_fn, str(i) + '.fna')
        yield open(fn, 'w'), i


def main():
    parser = argparse.ArgumentParser(
        description='Split a fasta file containing multiple sequences in multiple files containing one sequence'
    )
    parser.add_argument('output_dir', type=str, help='Path to the output directory')
    parser.add_argument(
        '-i', dest='input_file', metavar='STR', type=argparse.FileType('r'), default=sys.stdin,
        help='Fasta file [stdin]'
    )
    args = parser.parse_args()

    input_f = args.input_file
    output_dir_fn = args.output_dir

    start_seq = '>'
    split_f = split_fs(output_dir_fn)
    outfile = None

    for line in input_f:
        if start_seq not in line:
            outfile.write(line)
        else:
            if outfile:
                outfile.close()
            outfile, i = next(split_f)
            outfile.write(line)

    outfile.close()


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        # pipe error (e.g., when head is used)
        sys.stderr.close()
        exit(0)
    except KeyboardInterrupt:
        pro.message("Error: Keyboard interrupt")
        pro.close_log()
        exit(1)
