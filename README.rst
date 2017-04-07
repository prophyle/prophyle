ProPhyle – accurate and resource-frugal phylogeny-based metagenomic classification
==================================================================================


.. image:: https://travis-ci.org/karel-brinda/prophyle.svg?branch=master
	:target: https://travis-ci.org/karel-brinda/prophyle

ProPhyle is a metagenomic classifier based on BWT-index and phylogenetic trees,
whose indexing strategy is based on the bottom-up propagation of genomes' k-mers in the tree,
assembling contigs at each node and matching using a standard full-text search.
The analysis of shared k-mers between NGS reads and the genomes in the index determines
which nodes are the best candidates for their classification.


Getting started
---------------


Prerequisities
^^^^^^^^^^^^^^

* GCC 4.8+
* ZLib
* Python 3 with ete3 library
* SamTools


Installation using Conda (recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Environment installation::

	conda create -y --name prophyle -c etetoolkit -c bioconda \
		python==3.6 ete3 bitarray samtools=1.3.1
	source activate prophyle
	pip install --upgrade prophyle


Environment activation::

	source activate prophyle


Installation using PIP
^^^^^^^^^^^^^^^^^^^^^^

From PyPI::

	pip install --upgrade prophyle

From Git::

	pip install --upgrade git+https://github.com/karel-brinda/prophyle

From PyPI to the current directory::

	pip install --user prophyle
	export PYTHONUSERBASE=`pwd`
	export PATH=$PATH:`pwd`/bin


Pipeline example
^^^^^^^^^^^^^^^^

Quick example::

	prophyle download bacteria
	prophyle index -k 10 ~/prophyle/test_bacteria.nw test_idx
	prophyle classify test_idx reads.fq > result.sam
