.. _classify:


Classifying reads
=================

When the index construction is finished, you can classify your reads using

.. code-block:: bash

	prophyle classify <index_dir> <reads.fq>

For the index from above, the command would be:

.. code-block:: bash

	prophyle classify my_BV_index my_reads.fq

The assignments are reported in the SAM format.
Unless specified otherwise, the *k*-mer length is deduced from the index.
