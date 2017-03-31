.. _tutorials:

Tutorials
=========

Custom indexes
--------------

Building custom ProPhyle indexes is easy and fast. This example focuses on NCBI's data, but it is possible to build an index for virtually any set of sequences.

Download
^^^^^^^^

- Create a subdirectory in ProPhyle's main library directory, and change the working directory to the newly created one
- Download an ``assembly_summary.txt`` file from NCBI's FTP server. For RefSeq's bacterial genomes, follow this `link <ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt>`_. More information about these files `here <ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt>`_
- Select genomes of interest; the following command will select all complete genomes, but it is sufficient to add more conditions to ``awk`` to select a subset. Fields 1, 6 and 20 contain respectively the accession number, taxonomic identifier and ftp directory of each genome. ``acc2taxid`` will therefore be a tab separated file containing the taxid corresponding to each sequence's accession number. Please refer to `NCBI's genomes download FAQ <https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq>`_ for further information)::

    awk -F "\t" -v OFS="\t" '$12=="Complete Genome" && $11=="latest" {print $1, $6, $20}' assembly_summary.txt >ftpselection
    cut -f 3 ftpselection | awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' >ftpfilepaths
    cut -f 1,2 ftpselection >acc2taxid

- Download selected sequences using ``parallel --gnu -j 24 -a ftpfilepaths wget``. Number of jobs (``-j`` option) should be changed according to the number of available cores and bandwidth.
- Extract them using ``parallel --gnu -j 24 -a ftpfilepaths gunzip {/}``

Tree
^^^^

Build a taxonomic tree for the downloaded sequences using::

  ProPhyle_dir/bin/build_ncbi_taxtree.py -l <library_subdir> -t <acc2taxid_file> -o <output_file> -e <log_file>

Taxonomic identifiers are assigned to the sequences first, and then the tree is built using `ETE Toolkit <http://etetoolkit.org/>`_ and saved with newick format 1. Necessary node attributes are:

* ``name``: unique node name
* ``taxid``: unique taxonomic identifier
* ``seqname``: names of the sequences sharing the same taxid, separated by @
* ``fastapath``: paths of the sequences' fasta files, separated by @ (relative paths from their library directory)
* ``infasta_offset``: positions where each sequence starts inside the corresponding fasta files, separated by @
* ``base_len``: length of each sequence, separated by @

``seqname``, ``infasta_offset`` and ``base_len`` can be found in samtools' `fasta index <http://www.htslib.org/doc/faidx.html>`_.
Other optional attributes are ``sci_name``, ``named_lineage``, ``lineage``, ``rank`` (more info `here<http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html#automatic-tree-annotation-using-ncbi-taxonomy>`_). Using ETE library or modifying the above script it is possible to fastly adapt any tree to match the requirements above.

If taxonomic information is not available, ``taxid`` attribute can be replaced by the node name. This could be useful to perform classification over custom phylogenetic trees.

Index
^^^^^

Run the standard command to build ProPhyle index::

  prophyle index -n <tree.newick> -l <library_subdir> -k <kmer_length> <index_dir>
