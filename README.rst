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
can be found in our `poster <http://brinda.cz/publications/2017_prophyle_hitseq.pdf>`_.

Compared to other state-of-the-arts classifiers, ProPhyle provides several unique features:

* **Low memory requirements.** Compared to `Kraken <https://ccb.jhu.edu/software/kraken/>`_, ProPhyle has 9x smaller memory footprint for index construction and 5x smaller footprint for querying.
* **Flexibility.** ProPhyle is easy to use with any user-provided phylogenetic trees and reference genomes.
* **Standard bioinformatics formats.** `Newick/NHX <https://sites.google.com/site/cmzmasek/home/software/forester/nhx>`_ is used for representing phylogenetic trees and `SAM <https://samtools.github.io/hts-specs/SAMv1.pdf>`_ for reporting the assignments.
* **Lossless k-mer indexing.** ProPhyle stores a list of *all* genomes containing a *k*-mer.  It can be, therefore, accurate even with trees containing similar genomes (e.g, phylogenetic trees for a single species).
* **Deterministic behavior.** ProPhyle is a fully deterministic classifier with a mathematically well-defined behavior.

For information about how to use ProPhyle, see the main `ProPhyle documentation <http://prophyle.rtfd.io>`_.

Quick example
-------------

1. Clone the ProPhyle repository and add it to PATH: ::

        git clone --recursive http://github.com/karel-brinda/prophyle
        export PATH=$(pwd)/prophyle/prophyle:$PATH

2. Download the `RefSeq <https://www.ncbi.nlm.nih.gov/refseq/>`_ bacterial database: ::

        $ prophyle download bacteria

3. To quickly test ProPhyle functionality, create an index for randomly sampled 10% genomes from the E.coli subtree of NCBI taxonomy (with k=31): ::

        $ prophyle index -s 0.1 ~/prophyle/bacteria.nw@561 _index_ecoli

4. Classify your reads: ::

        $ prophyle classify _index_ecoli reads.fq > result.sam
