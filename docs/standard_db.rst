.. _standard_db:

Building standard databases
===========================

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





