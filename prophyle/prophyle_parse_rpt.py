#! /usr/bin/env python3

import os
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="")
    parser.add_argument(type=str, dest='library_dir')

    args = parser.parse_args()

    for dirpath, _, filenames in os.walk(args.library_dir):
        for filename in (f for f in filenames if f.endswith('.rpt')):
            fn = os.path.join(dirpath, filename)
            with open(fn, 'r') as in_rpt:
                for line in in_rpt:
                    if line.startswith('Accession:'):
                        acc = (line.split(':')[1]).split('.')[0].strip()
                    elif line.startswith('Taxid:'):
                        tax = line.split(':')[1].strip()
                        print(acc, tax, sep='\t')
