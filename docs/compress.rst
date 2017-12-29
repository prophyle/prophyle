.. _compress:

Compressing indexes for transmission
====================================

To simplify archiving indexes and their transmission over the Internet, ProPhyle
provides a compression feature. This compression is based
on creating a ``.tar.gz`` archive containing all information necessary to
reconstruct the entire index.
For a specification of the format, see Section :doc:`formats`.


Compression
-----------

An index can be compressed by

.. code-block:: bash

	prophyle compress <index.dir> [<archive.tar.gz>]

If ``<archive.tar.gz>`` is not provided, ProPhyle will derive the filename from the name of the index.



Decompression
-------------

A compressed index can then be decompressed by

.. code-block:: bash

	prophyle decompress [-K] <archive.tar.gz> [<output.dir>]

Decompression includes reconstructing the sampled suffix array and the k-LCP bit array.
If the option ``-K`` is provided, k-LCP will not be reconstructed.
