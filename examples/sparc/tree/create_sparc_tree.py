#! /usr/bin/env python3

import urllib.request
import os
import sys
import ete3

PUBLISHED_TREE="http://datadryad.org/bitstream/handle/10255/dryad.83423/SPARC.core_genes.tree"
OUTPUT_NW="sparc.nw"

def create_sparc_tree(tree_url,tree_fn):

	nw_source= urllib.request.urlopen(tree_url).read().decode('utf-8')
	#print(nw_source.decode('utf-8'))

	tree=ete3.Tree(nw_source,format=1)

	i=1
	for node in tree.traverse("postorder"):
		if node.name=="":
			nname="int_{}".format(i)
			node.name=nname
			i+=1
		else:
			node.add_features(fastapath="{}.fa".format(node.name))

	tree.write(
			format=1,
			features=['fastapath'],
			outfile=tree_fn,
			format_root_node=True,
		)

def main():
	create_sparc_tree(
		PUBLISHED_TREE,
		OUTPUT_NW,
	)


if __name__ == "__main__":
	main()
