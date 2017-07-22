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

* `kraken-report <https://ccb.jhu.edu/software/kraken/MANUAL.html#sample-reports>`_ format
* `MetaPhlAn2 <https://bitbucket.org/biobakery/biobakery/wiki/metaphlan2#rst-header-output-files>`_ format
* Custom `Centrifuge <https://ccb.jhu.edu/software/centrifuge/manual.shtml#centrifuge-summary-output-the-default-filename-is-centrifuge_report.tsv>`_ format

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
