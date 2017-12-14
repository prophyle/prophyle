:tocdepth: 1

.. _classification:

Read classification
===================

When the index construction is finished, you can classify your reads using ::

        prophyle classify <index_dir> <reads.fq>

For the index from above, the command would be: ::

        prophyle classify my_BV_index my_reads.fq

The assignments are reported in the SAM format.
Unless specified otherwise, the *k*-mer length is deduced from the index.
