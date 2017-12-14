:tocdepth: 1

.. _example:

Quick example
=============

Examples
--------

Quick test (small k, subsampled bacterial database)::

	prophyle download bacteria
	prophyle index -k 10 -s 0.1 ~/prophyle/bacteria.nw test_idx
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
