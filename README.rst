ProPhyle - accurate and resource-frugal phylogeny-based metagenomic classification
==================================================================================


.. image:: https://travis-ci.org/karel-brinda/prophyle.svg?branch=master
	:target: https://travis-ci.org/karel-brinda/prophyle

ProPhyle is a metagenomic classifier based on BWT-index and phylogenetic trees, whose indexing strategy is based on the bottom-up propagation of genomes' k-mers in the tree, assembling contigs at each node and matching using a standard full-text search. The analysis of shared k-mers between NGS reads and the genomes in the index determines which nodes are the best candidates for their classification.

Getting started
---------------

Prerequisities
^^^^^^^^^^^^^^

* GCC 4.8+
* ZLib
* Python 3 with ete3 library
* SamTools

Recommended way of installation using Conda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Environment installation::

	conda create -y --name prophyle \
		-c etetoolkit -c bioconda \
		python==3.4 ete3 bitarray \
		parallel blast samtools=1.3.1


Environment activation::

	source activate prophyle


Compile all programs::

  make

Pipeline example
^^^^^^^^^^^^^^^^

prophyle download bacteria
prophyle index -t bacteria.nw -g bacteria/ idx_bacteria
prophyle classify idx_bacteria reads.fq > result.sam
