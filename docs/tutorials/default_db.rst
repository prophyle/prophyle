.. _default_db:

Default databases
=================


Download
--------

Download libraries using ``prophyle download <library>``. Use ``all`` to
download bacteria, viruses and plasmids from NCBI's RefSeq archive and
`HMP <http://hmpdacc.org/>`_ sequences.

You can find pre-built trees in the ``trees`` directory (they are automatically
copied to ProPhyle's home directory when a library is downloaded).


Index
-----

Run ``prophyle index -k <kmer_length> <tree_1.nw> [<tree_2.nw> ...] <index_dir>``
to build an index for a subsection of the above libraries. If multiple trees are
given they are merged automatically using the TaxID information of their
nodes.


Classify
--------

Run ``prophyle classify <index_dir> <reads.fq>`` to classify reads in ``fastq``
format and output the results in SAM format to the standard output.


Examples
--------

Quick test (small k, subsampled bacterial database)::

	prophyle download bacteria
	prophyle index -k 10 ~/prophyle/test_bacteria.nw test_idx
	prophyle classify test_idx reads.fq > result.sam

Bacterial database (k=31)::

	prophyle download bacteria
	prophyle index -k 31 ~/prophyle/bacteria.nw idx_bac
	prophyle classify idx_bac reads.fq > result.sam

Example for bacterial and viral database (k=31)::

	prophyle download bacteria
	prophyle download viruses
	prophyle index -k 31 ~/prophyle/bacteria.nw ~/prophyle/viruses.nw idx_bac_vir
	prophyle classify idx_bac_vir reads.fq > result.sam
