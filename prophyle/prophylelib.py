#! /usr/bin/env python3

"""Auxiliary ProPhyle functions.

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT
"""

import sys
import ete3

def load_nhx_tree(nhx_fn, validate=True):
	tree=ete3.Tree(
		nhx_fn,
		format=1
	)

	if validate:
		validate_prophyle_nhx_tree(tree)

	return tree


def save_nhx_tree(tree, nhx_fn):
	assert isinstance(tree, ete3.Tree)

	# make saving newick reproducible
	features=set()
	for n in tree.traverse():
		features|=n.features

	# otherwise some names stored twice – also as a special attribute
	features.remove("name")

	tree.write(
		features=sorted(features),
		format = 1,
		format_root_node = True,
		outfile = nhx_fn,
	)


def validate_prophyle_nhx_tree(tree, verbose=True, throw_exceptions=True, output_file=sys.stderr):

	node_names=set()

	error=False

	existing_names=[]
	underscored_names=[]

	without_name=[]
	empty_name=[]

	duplicates=[]

	for i,node in enumerate(tree.traverse("postorder")):
		noname=True
		try:
			nname=node.name
			noname=False
		except AttributeError:
			without_name.append((i,None))

		if not noname:
			if nname=='':
				empty_name.append((i,nname))
				error=True
			if nname in existing_names:
				duplicates.append((i,nname))
				error=True
			if "@" in nname:
				underscored_names.append((i,nname))
				error=True
			existing_names.append((i,nname))

	def _format_node_list(node_list):
		return ", ".join(map(str,node_list))

	def _error_report(node_list, message):
		if len(node_list)>0:
			print("   * {} nodes {}: {}".format(
					len(node_list),
					message,
					_format_node_list(node_list),
				),file=output_file)

	if verbose:
		if error:
			print("Error",file=output_file)

		_error_report(without_name, "without name")
		_error_report(empty_name, "with empty name")
		_error_report(duplicates, "with a duplicate name")
		_error_report(underscored_names, "a name containing '_'")

	if throw_exceptions:
		if error:
			raise ValueError("Invalid tree. The format of the tree is not correct. See the messages above.")

	if error:
		return False
	else:
		return True
