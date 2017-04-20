ProPhyle – accurate and resource-frugal phylogeny-based metagenomic classification
==================================================================================


.. image:: https://travis-ci.org/karel-brinda/prophyle.svg?branch=master
	:target: https://travis-ci.org/karel-brinda/prophyle

.. image:: https://readthedocs.org/projects/prophyle/badge/?version=latest
	:target: http://prophyle.rtfd.org
	
.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square
	:target: https://anaconda.org/bioconda/prophyle

.. image:: https://badge.fury.io/py/prophyle.svg
    :target: https://badge.fury.io/py/prophyle

ProPhyle is a metagenomic classifier based on BWT-index and phylogenetic trees.
Its indexing strategy uses a bottom-up propagation of k-mers in the tree,
assembling contigs at each node, and matching using a standard full-text search.
The analysis of shared k-mers between NGS reads and the genomes in the index determines
which nodes are the best candidates for their classification.

More information can be found in our `poster <http://brinda.cz/publications/2017_prophyle_hsph_poster_day.pdf>`_.


Getting started
---------------

Prerequisities
^^^^^^^^^^^^^^

* GCC 4.8+
* ZLib
* Python 3 with ete3 library
* SamTools



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


Installation using Bioconda
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Environment installation::

	conda create -c bioconda -n prophyle prophyle

Environment activation::

	source activate prophyle


Examples
^^^^^^^^

Quick test (small k, subsampled bacterial database)::

	prophyle download bacteria
	prophyle index -k 10 ~/prophyle/test_bacteria.nw test_idx
	prophyle classify test_idx reads.fq > result.sam

Bacterial database (k=31)::

	prophyle download bacteria
	prophyle index -k 31 ~/prophyle/bacteria.nw idx_bac
	prophyle classify idx_bac reads.fq > result.sam

Bacterial and viral database (k=31)::

	prophyle download bacteria viruses
	prophyle index -k 31 ~/prophyle/bacteria.nw ~/prophyle/viruses.nw idx_bac_vir
	prophyle classify idx_bac_vir reads.fq > result.sam
