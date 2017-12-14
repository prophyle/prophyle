:tocdepth: 1

ProPhyle
========

ProPhyle is a metagenomic classifier based on BWT-index and phylogenetic trees,
whose indexing strategy is based on the bottom-up propagation of genomes’
k-mers in the tree, assembling contigs at each node and matching using a
standard full-text search. The analysis of shared k-mers between NGS reads and
the genomes in the index determines which nodes are the best candidates for
their classification.

Links
-----

`GitHub repository`_ -
`Bug reporting`_ -
`Contact`_

.. _`GitHub repository`: http://github.com/karel-brinda/ProPhyle
.. _`Bug reporting`: http://github.com/karel-brinda/ProPhyle/issues
.. _`Contact`: 	kbrinda@hsph.harvard.edu


Table of content
----------------


.. toctree::
	:maxdepth: 3
	:name: mastertoc

	index
	requirements
	install
	databases
	classification
	formats
	reference
