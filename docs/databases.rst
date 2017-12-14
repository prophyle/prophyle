:tocdepth: 1

.. _databases:

Building databases
==================

ProPhyle comes with several genome libraries containing
RefSeq genomes, augmented with the NCBI taxonomy.

Downloading genomes
-------------------

These libraries can be downloaded using ``prophyle download <library> [<library> ...]``,
where ``<library>`` should be replaced by ``bacteria``, ``viruses``, or ``plasmids``.
The command also copies a prebuild Newick/NHX tree for the specified library.
If the `-d` parameter is not specified, all files are placed to `~/prophyle`.


To download all viral and bacterial genomes from RefSeq, execute: ::

        prophyle download bacteria viruses


.. Use ``all`` to

.. download bacteria, viruses and plasmids from NCBI's RefSeq archive and

..  `HMP <http://hmpdacc.org/>`_ sequences.

..

.. You can find pre-built trees in the ``trees`` directory (they are automatically

.. copied to ProPhyle's home directory when a library is downloaded).


Index construction
------------------

Once a library is downloaded, a ProPhyle index can be constructed using: ::

        prophyle index [-g DIR] [-j INT] [-k INT] [-M] [-P] [-K] <tree.nw> [<tree.nw> ...] <index.dir>

`<tree.nw>` is a Newick/NHX for the index. The trees from the previous command
are placed in `~/prophyle` and they are called `bacteria.nw`, `viruses.nw`, etc.
`<index.dir>` is the directory directory where your index files are going to
be placed.

There are multiple other parameters that can be used.
`-j` can be used to specify the number of CPU cores used for index construction (all cores are used otherwise).
`-k` serves to set the *k*-mer length (31 in default).
`-M` activates low complexity regions filtering using DustMasker. Please, ensure that the program is install (try to run `dustmasker`).
If multiple trees are used, they are going to be merged. Therefore, a name collision can
appear. To prevent such a situation, ProPhyle prepends numerical prefixes to the
node names (unless `-P` is used).
The `-K` parameter can be used to deactivate *k*-LCP array construction. The resulting index
would be slightly smaller, but querying would become much slower.

So the entire command for index construction can look, for instance,
like this: ::

        prophyle index -k 25 ~/prophyle/bacteria.nw ~/prophyle/viruses.nw my_BV_index

Index construction might take several hours, based on the database size, *k* and the number
of used cores.


Classification of reads
-----------------------

When the index construction is finished, you can classify your reads using ::

        prophyle classify <index_dir> <reads.fq>

For the index from above, the command would be: ::

        prophyle classify my_BV_index my_reads.fq

The assignments are reported in the SAM format.
Unless specified otherwise, the *k*-mer length is deduced from the index.


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
