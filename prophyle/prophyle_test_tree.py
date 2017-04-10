#! /usr/bin/env python3

import os
import sys
import argparse

from ete3 import Tree

import logging

DEFAULT_FORMAT = 1


def verify_newick_tree(newick):

	node_names=set()

	error=False

	existing_names=set()
	underscored_names=set()
	without_name=0
	empty_name=0
	duplicates=[]

	for node in tree.traverse("postorder"):
		noname=True
		try:
			rname=node.name
			noname=False
		except AttributeError:
			without_name+=1

		if not noname:
			if rname=='':
				empty_name+=1
				error=True
			if rname in existing_names:
				duplicates.append(rname)
				error=True
			if "_" in rname:
				underscored_names.add(underscored_names)
				error=True
			existing_names.add(rname)


	duplicates.sort()

	if error:
		print("Error",file=sys.stderr)

	if without_name>0:
		print("{} nodes without name".format(without_name),file=sys.stderr)

	if empty_name>0:
		print("{} nodes with empty name".format(empty_name),file=sys.stderr)

	if len(duplicates)>0:
		print("{} node(s) with a duplicate name: {}".format(len(duplicates), ", ".join(duplicates)),file=sys.stderr)

	if len(underscored_names)>0:
		print("{} node(s) with a name containing '_': {}".format(len(underscored_names), ", ".join(underscored_names)),file=sys.stderr)

	if error:
		sys.exit(1)


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Verify a Newick tree')

	parser.add_argument('-n', '--newick-tree',
			type=str,
			metavar='str',
			required=True,
			dest='newick_fn',
			help='newick tree',
		)

	parser.add_argument('-p', '--print-tree',
			action='store_true',
			dest='printt',
			help='print the newick tree',
		)

	args = parser.parse_args()
	newick_fn=args.newick_fn
	p=args.printt



	tree=Tree(newick_fn,format=DEFAULT_FORMAT)

	if p:
		print(tree.get_ascii(show_internal=True))

	verify_newick_tree(tree)
