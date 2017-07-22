.. _reference:

Reference
=========

Here you can find format specifications and an automatically generated reference
for ProPhyle's CLI.

Formats
-------

Trees
^^^^^

`Newick <http://evolution.genetics.washington.edu/phylip/newicktree.html>`_
format `1 <http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-and-writing-newick-trees>`_
with NHX annotations, which can be easily created and modified using the
`ete3 <http://etetoolkit.org/>`_ python package.

Classification output
^^^^^^^^^^^^^^^^^^^^^

Support for both `SAM <http://samtools.github.io/hts-specs/>`_ and `Kraken <https://ccb.jhu.edu/software/kraken/MANUAL.html#output-format>`_
output formats.

Analysis output
^^^^^^^^^^^^^^^

* `kraken-report <https://ccb.jhu.edu/software/kraken/MANUAL.html#sample-reports>`_ format:
	1. Percentage of reads covered by the clade rooted at this taxon
	2. Number of reads covered by the clade rooted at this taxon
	3. Number of reads assigned directly to this taxon
	4. A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are simply '-'.
	5. NCBI taxonomy ID
	6. indented scientific name
* `MetaPhlAn2 <https://bitbucket.org/biobakery/biobakery/wiki/metaphlan2#rst-header-output-files>`_ format:
	1. clades, ranging from taxonomic kingdoms (Bacteria, Archaea, etc.) through species.
	The taxonomic level of each clade is prefixed to indicate its level: Kingdom: ``k__``, Phylum: ``p__``, Class: ``c__``, Order: ``o__``, Family: ``f__``, Genus: ``g__``, Species: ``s__``.
	Since sequence-based profiling is relative and does not provide absolute cellular abundance measures,
	clades are hierarchically summed. Each level will sum to 100%; that is, the sum of all kindom-level
	clades is 100%, the sum of all genus-level clades (including unclassified) is also 100%, and so forth.
	OTU equivalents can be extracted by using only the species-level ``s__`` clades from this file
	(again, making sure to include clades unclassified at this level).
* Custom `Centrifuge <https://ccb.jhu.edu/software/centrifuge/manual.shtml#centrifuge-summary-output-the-default-filename-is-centrifuge_report.tsv>`_ format::

	#name                                                           taxID   taxRank    kmerCount   numReads   numUniqueReads   abundance
	Wigglesworthia glossinidia endosymbiont of Glossina brevipalpis 36870   leaf       703004      5981.37    5964             0

1. name of a genome, or the name corresponding to a taxonomic ID (the second column) at a rank higher than the strain (e.g., Wigglesworthia glossinidia endosymbiont of Glossina brevipalpis).
2. taxonomic ID (e.g., 36870).
3. taxonomic rank (e.g., leaf).
4. number of k-mers propagated up to the node (e.g., 703004).
5. number of reads classified to this node including multi-classified reads (divided by the number of assignments, e.g., 5981.37).
6. number of reads uniquely classified to this genomic sequence (e.g., 5964).
7. not used yet.

Main program's reference
------------------------

``prophyle`` (list of subcommands)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	.. include:: main.txt
		:code: text

``prophyle download``
^^^^^^^^^^^^^^^^^^^^^
	.. include:: download.txt
		:code: text

``prophyle index``
^^^^^^^^^^^^^^^^^^
	.. include:: index.txt
		:code: text

``prophyle classify``
^^^^^^^^^^^^^^^^^^^^^
	.. include:: classify.txt
		:code: text

``prophyle analyze``
^^^^^^^^^^^^^^^^^^^^
	.. include:: analyze.txt
		:code: text

``prophyle compress``
^^^^^^^^^^^^^^^^^^^^^
	.. include:: compress.txt
		:code: text

``prophyle decompress``
^^^^^^^^^^^^^^^^^^^^^^^
	.. include:: decompress.txt
		:code: text

Other programs' reference
-------------------------

``prophyle_ncbi_tree``
^^^^^^^^^^^^^^^^^^^^^^
	.. include:: prophyle_ncbi_tree.txt
		:code: text

``prophyle_assembler``
^^^^^^^^^^^^^^^^^^^^^^
	.. include:: prophyle_assembler.txt
		:code: text

``prophyle_index`` (list of subcommands)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	.. include:: prophyle_index.txt
		:code: text

``prophyle_index build``
^^^^^^^^^^^^^^^^^^^^^^^^
.. include:: prophyle_index_build.txt
	:code: text

``prophyle_index query``
^^^^^^^^^^^^^^^^^^^^^^^^
.. include:: prophyle_index_query.txt
	:code: text

``prophyle_assignment``
^^^^^^^^^^^^^^^^^^^^^^^
	.. include:: prophyle_assignment.txt
		:code: text

``prophyle_analyze``
^^^^^^^^^^^^^^^^^^^^
	.. include:: prophyle_analyze.txt
		:code: text
