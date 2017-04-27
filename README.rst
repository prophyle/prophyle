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


Introduction
------------

ProPhyle is a *k*-mer based metagenomic classifier using Burrows-Wheeler Transform.
Its indexing strategy relies on a bottom-up propagation of *k*-mers in the tree,
assembling contigs at each node, and matching using a standard full-text search using BWT-index.
The analysis of shared *k*-mers between NGS reads and the genomes in the index determines
which nodes are the best candidates for their classification.
More information about the indexing scheme
can be found in our `poster <http://brinda.cz/publications/2017_prophyle_hsph_poster_day.pdf>`_.

Compared to other state-of-the-arts classifiers, ProPhyle provides several unique features:

* **Low memory requirements.** Compared to Kraken, ProPhyle has 9x smaller memory footprint for index construction and 5x smaller footprint for querying.
* **Flexibility.** ProPhyle is easy to use with any user-provided phylogenetic tree or multiple trees.
* **Standard bioinformatics formats.** Newick/NHX is used for representing phylogenetic trees and SAM for reporting the assignments.
* **Lossless k-mer indexing.** ProPhyle stores a list of *all* genomes containing a *k*-mer.
  It can be, therefore, accurate even with trees containing similar genomes
  (e.g, phylogenetic trees for a single species).

For information about how to use ProPhyle, see the main `ProPhyle documentation <http://prophyle.rtfd.io>`_.

Quick example
-------------

1. Create a Bioconda environment with ProPhyle and activate it: ::

        $ conda create -c bioconda -n prophyle prophyle
        $ source activate prophyle

2. To quickly test ProPhyle functionality, download the bacterial database, and create an index
with *k*-mer legth 12 and a subset of bacterial genomes: ::

        $ prophyle index -k 10 ~/prophyle/test_bacteria.nw test_idx
	$ prophyle download bacteria

3. Classify your reads: ::

	prophyle classify test_idx reads.fq > result.sam

