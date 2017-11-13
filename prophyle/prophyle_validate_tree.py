#! /usr/bin/env python3
"""Test whether given Newick/NHX trees are valid for ProPhyle.

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT

Example:

    $ prophyle_validate_tree.py ~/prophyle/bacteria.nw ~/prophyle/viruses.nw
"""

import os
import sys
import argparse

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro


def main():
    parser = argparse.ArgumentParser(description='Verify a Newick/NHX tree')

    parser.add_argument(
        'tree',
        metavar='<tree.nw>',
        type=str,
        nargs='+',
        help='phylogenetic tree (in Newick/NHX)',
    )

    args = parser.parse_args()
    tree_fns = args.tree

    ok = True

    for tree_fn in tree_fns:
        print("Validating '{}'".format(tree_fn))
        tree = pro.load_nhx_tree(tree_fn, validate=False)
        r = pro.validate_prophyle_nhx_tree(tree, verbose=True, throw_exceptions=False, output_fo=sys.stdout)
        if r:
            print("   ...OK")
        else:
            ok = False
        print()

    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
