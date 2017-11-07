#! /usr/bin/env python3
"""Plot a Newick/NHX tree.

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT

Example:

    Plot tree:

        $ prophyle_plot_tree.py ~/prophyle/bacteria.nw tree.pdf

    Plot tree of scientific names:

        $ ./prophyle_plot_tree.py -a sci_name ~/prophyle/bacteria.nw sci_names.pdf
"""

import os
import sys
import argparse
import ete3
import locale

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro


def plot_tree(newick_in_fn, out_plot_fn, attribute_name):
    tree = pro.load_nhx_tree(newick_in_fn)

    ts = ete3.TreeStyle()
    ts.show_leaf_name = False

    def my_layout(node):
        name = getattr(node, attribute_name)

        try:
            kmer_full = locale.format("%d", int(node.kmers_full), grouping=True),
        except AttributeError:
            kmer_full = None

        try:
            kmer_reduced = locale.format("%d", int(node.kmers_reduced), grouping=True)
        except AttributeError:
            kmer_reduced = None

        if kmer_full is None:
            if kmer_reduced is None:
                t = name
            else:
                t = "{} [red. {}]".format(name, kmer_reduced)
        else:
            if kmer_reduced is None:
                t = "{} [full {}]".format(name, kmer_full)
            else:
                t = "{} [full {} & red. {}]".format(name, kmer_full, kmer_reduced)

        f = ete3.TextFace(t, tight_text=True)
        ete3.add_face_to_node(f, node, column=0, position="branch-right")

    ts.layout_fn = my_layout
    tree.render(out_plot_fn, tree_style=ts)


def main():
    parser = argparse.ArgumentParser(description='Plot a Newick/NHX tree')

    parser.add_argument(
        'newick_in_fn',
        metavar='<tree.nhx>',
        type=str,
        help='phylogenetic tree (in Newick/NHX)',
    )

    parser.add_argument(
        'out_plot_fn',
        metavar='<figure.{pdf,png,svg,..}>',
        type=str,
        help='output figure',
    )

    parser.add_argument(
        '-a',
        metavar='str',
        dest='attribute_name',
        type=str,
        help='attribute to print with each node (e.g., sci_name)',
        default='name',
    )

    args = parser.parse_args()

    plot_tree(
        newick_in_fn=args.newick_in_fn,
        out_plot_fn=args.out_plot_fn,
        attribute_name=args.attribute_name,
    )


if __name__ == "__main__":
    main()
