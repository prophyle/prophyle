ProPhyle - accurate and resource-frugal phylogeny-based metagenomic classification
==================================================================================


.. image:: https://travis-ci.org/karel-brinda/prophyle.svg?branch=master
	:target: https://travis-ci.org/karel-brinda/prophyle

Getting started
---------------

Prerequisities
~~~~~~~~~~~~~~

* GCC 4.8+
* ZLib
* Python 3 with ete3 library
* SamTools

Recommended way of installation using Conda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Environment installation:

.. code-block:: bash

	conda create -y --name prophyle \
		-c etetoolkit -c bioconda \
		python==3.4 ete3 bitarray \
		parallel blast samtools=1.3.1


Environment activation:

.. code-block:: bash

        source activate prophyle


Compile all programs

.. code-block:: bash

  make -C src

Custom taxonomic trees
~~~~~~~~~~~~~~~~~~~~~~

Use `bin/build_taxonomic_tree.py <bin/build_taxonomic_tree.py>`_ to build custom taxonomic trees starting from your database's fasta indexes and taxonomy files (`library/Taxonomy <library/Taxonomy>`_ for more information). Taxonomic identifiers are assigned to the sequences first, and then the tree is built using `ETE Toolkit <http://etetoolkit.org/>`_ and saved with newick format 1. Necessary node attributes are:

* ``name``: unique node name
* ``taxid``: unique taxonomic identifier
* ``seqname``: names of the sequences sharing the same taxid, separated by @
* ``fastapath``: paths of the sequences' fasta files, separated by @ (relative paths from their library directory)
* ``infasta_offset``: positions where each sequence starts inside the corresponding fasta files, separated by @
* ``base_len``: length of each sequence, separated by @

``seqname``, ``infasta_offset`` and ``base_len`` can be found in samtools' `fasta index <http://www.htslib.org/doc/faidx.html>`_.
Other optional attributes are ``sci_name``, ``named_lineage``, ``lineage``, ``rank`` (more info `here <http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html#automatic-tree-annotation-using-ncbi-taxonomy>`_).