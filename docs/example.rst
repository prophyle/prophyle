:tocdepth: 1

.. _example:

Quick example
=============

1. Install ProPhyle using `Bioconda <https://bioconda.github.io/>`_:

	.. code-block:: bash

		conda config --add channels defaults
		conda config --add channels conda-forge
		conda config --add channels bioconda
		conda install prophyle

2. Download the `RefSeq <https://www.ncbi.nlm.nih.gov/refseq/>`_ bacterial database:

	.. code-block:: bash

		prophyle download bacteria

3. To quickly test ProPhyle functionality, create an index for randomly sampled 10% genomes from the E.coli subtree of the NCBI taxonomy with k-mer size 25:

	.. code-block:: bash

		prophyle index -s 0.1 -k 25 ~/prophyle/bacteria.nw@561 _index_ecoli

4. Classify your reads:

	.. code-block:: bash

		prophyle classify _index_ecoli reads.fq > result.sam

