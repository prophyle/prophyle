.. _formats:


File formats
============

.. contents::
	:depth: 3
	:local:
	:backlinks: none


Input data formats
------------------

Newick trees
^^^^^^^^^^^^

`Newick <http://evolution.genetics.washington.edu/phylip/newicktree.html>`_ trees are eventually transformed to NHX.

Examples
""""""""

	A Newick tree with named leaves::

		((n1,n2,n3),(n5,n6));


	A Newick tree with nameds nodes::

		((n1,n2)o1,(n3,n4,n5)o2)p1;



NHX trees
^^^^^^^^^

`New Hampshire X Format <https://sites.google.com/site/cmzmasek/home/software/forester/nhx>`_
is parsed using the `ETE3 library <http://etetoolkit.org/>`_  (see specification of `format 1 <http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-and-writing-newick-trees>`_).


	.. list-table:: NHX attributes
	   :widths: 7 7 20
	   :header-rows: 1

	   * - Attribute
	     - Type
	     - Description
	   * - name
	     - string
	     - unique node name (ideally the TaxID of the node)
	   * - taxid
	     -
	     - unique taxonomic identifier
	   * - seqname
	     -
	     - names of the sequences sharing the same taxid, separated by @
	   * - fastapath
	     -
	     - paths of the sequences' fasta files, separated by @ (relative paths from ProPhyle's home directory)
	   * - base_len
	     -
	     - length of each sequence, separated by @


Example
"""""""
	Previous tree after autocomplete::

		(((n1:1[&&NHX:dist=1.0:fastapath=n1.fa:support=1.0],n2:1[&&NHX:dist=1.0:fastapath=n2.fa:support=1.0])o1:1[&&NHX:dist=1.0:support=1.0],(n3:1[&&NHX:dist=1.0:fastapath=n3.fa:support=1.0],n4:1[&&NHX:dist=1.0:fastapath=n4.fa:support=1.0],n5:1[&&NHX:dist=1.0:fastapath=n5.fa:support=1.0])o2:1[&&NHX:dist=1.0:support=1.0])p1:0[&&NHX:dist=0.0:support=1.0])merge_root:1[&&NHX:dist=1.0:support=1.0];


Sequences
^^^^^^^^^

Input sequences can be provided in the FASTA or FASTQ formats. Any non-ACGT characters are treated as
unknown nucleotides. All k-mers containing an unknown nucleotide are discarded.


.. <hr />

Assignments
-----------

Read assignments in SAM/BAM
^^^^^^^^^^^^^^^^^^^^^^^^^^^

`SAM format <http://samtools.github.io/hts-specs/>`_


	.. list-table:: SAM fields
	   :widths: 3 3 20
	   :header-rows: 1

	   * - Column
	     - Name
	     - Description
	   * - 1
	     - QNAME
	     - Query name
	   * - 2
	     - FLAG
	     - ``0`` if assigned, ``4`` otherwise
	   * - 3
	     - RNAME
	     - Node name
	   * - 4
	     - POS
	     - ``1`` if assigned, unused (``0``) otherwise
	   * - 5
	     - MAPQ
	     - ``60`` if assigned, unused (``0``) otherwise
	   * - 6
	     - CIGAR
	     - Coverage bit-mask (e.g., `7=3X3=` means `1111111000111`) if assigned, unused (``*``) otherwise
	   * - 7
	     - RNEXT
	     - Unused (``*``)
	   * - 8
	     - PNEXT
	     - Unused (``0``)
	   * - 9
	     - TLEN
	     - Unused (``0``)
	   * - 10
	     - SEQ
	     - Sequence of bases if ``-P``, unused (``*``) otherwise
	   * - 11
	     - QUAL
	     - Base qualities if ``-P``, unused (``*``) otherwise


	.. list-table:: SAM tags
	   :widths: 5 5 20
	   :header-rows: 1

	   * - Tag
	     - Type
	     - Description
	   * - h1
	     - integer
	     - Number of shared k-mers
	   * - h2
	     - float
	     - Normalized number of shared k-mers
	   * - hf
	     - float
	     - Proportion of hits
	   * - c1
	     - integer
	     - Number of covered positions in the query
	   * - c2
	     - float
	     - Normalized number of covered positions in the query
	   * - cf
	     - float
	     - Proportion of coverage
	   * - is
	     - int
	     - Number of reported assignments for the query
	   * - ii
	     - int
	     - ID of the curent assignment
	   * - hc
	     - string
	     - Hit bit-mask  (e.g., `7=1X3=` means `11111110111`)


Read assignments in a Kraken-like format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Kraken <https://ccb.jhu.edu/software/kraken/MANUAL.html#output-format>`_


Abundances estimates (experimental)
-----------------------------------

Abundances in the Kraken report format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`kraken-report <https://ccb.jhu.edu/software/kraken/MANUAL.html#sample-reports>`_ format:


	.. list-table:: Kraken report format
	   :widths: 5 20
	   :header-rows: 1

	   * - Column
	     - Description
	   * - 1
	     - Percentage of reads covered by the clade rooted at this taxon
	   * - 2
	     - Number of reads covered by the clade rooted at this taxon
	   * - 3
	     - Number of reads assigned directly to this taxon
	   * - 4
	     - A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are simply '-'.
	   * - 5
	     - NCBI taxonomy ID
	   * - 6
	     - Indented scientific name


Abundances in the Metaphlan2 report format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`MetaPhlAn2 <https://bitbucket.org/biobakery/biobakery/wiki/metaphlan2#rst-header-output-files>`_ format

	.. list-table:: Metaphlan 2 report format
	   :widths: 5 20
	   :header-rows: 1

	   * - Column
	     - Description
	   * - 1
	     - Clades, ranging from taxonomic kingdoms (Bacteria, Archaea, etc.) through species
	   * - 2
	     - The taxonomic level of each clade is prefixed to indicate its level: Kingdom: ``k__``, Phylum: ``p__``, Class: ``c__``, Order: ``o__``, Family: ``f__``, Genus: ``g__``, Species: ``s__``



Since sequence-based profiling is relative and does not provide absolute cellular abundance measures, clades are hierarchically summed. Each level will sum to 100%; that is, the sum of all kindom-level clades is 100%, the sum of all genus-level clades (including unclassified) is also 100%, and so forth. OTU equivalents can be extracted by using only the species-level ``s__`` clades from this file (again, making sure to include clades unclassified at this level).


Abundances in the Centrifuge report format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Centrifuge <https://ccb.jhu.edu/software/centrifuge/manual.shtml#centrifuge-summary-output-the-default-filename-is-centrifuge_report.tsv>`_ format.

Specification
"""""""""""""

	.. list-table:: Centrifuge format
	   :widths: 5 20
	   :header-rows: 1


	   * - Column
	     - Description
	   * - 1
	     - name of a genome, or the name corresponding to a taxonomic ID (the second column) at a rank higher than the strain (e.g., Wigglesworthia glossinidia endosymbiont of Glossina brevipalpis).
	   * - 2
	     - taxonomic ID (e.g., 36870).
	   * - 3
	     - taxonomic rank (e.g., leaf).
	   * - 4
	     - number of k-mers propagated up to the node (e.g., 703004).
	   * - 5
	     - number of reads classified to this node including multi-classified reads (divided by the number of assignments, e.g., 5981.37)
	   * - 6
	     - number of reads uniquely classified to this genomic sequence (e.g., 5964)
	   * - 7
	     - unused


Example
"""""""

::

	#name                                                           taxID   taxRank    kmerCount   numReads   numUniqueReads   abundance
	Wigglesworthia glossinidia endosymbiont of Glossina brevipalpis 36870   leaf       703004      5981.37    5964             0


Internal ProPhyle formats
-------------------------

ProPhyle Index
^^^^^^^^^^^^^^

Directory with the following files:

	.. list-table:: ProPhyle index
	   :widths: 5 20
	   :header-rows: 1

	   * - File name
	     - Description
	   * - ``index.fa``
	     - Assembled contigs, name of sequences are of the following format: ``[node_name]@c[contig_id]``
	   * - ``index.fa.amb``
	     - List of ambiguous nucleotides, no values
	   * - ``index.fa.ann``
	     - List of contigs and their starting positions in the master string
	   * - ``index.fa.[k].klcp``
	     - k-LCP array
	   * - ``index.fa.bwt``
	     - Burrows-Wheeler Transform of the master string (merged sequences + reverse completement) + OCC table (BWA format)
	   * - ``index.fa.kmers.tsv``
	     - k-mer statistics, format: ``[node_name].[full|reduced].fa	[#kmers]``, where ``full`` refers to all associated k-mers and ``reduced`` to represented k-mers
	   * - ``index.fa.pac``
	     - Packed sequences (BWA format)
	   * - ``index.fa.sa``
	     - Sampled suffix array (BWA format)
	   * - ``index.json``
	     - Index parameters: k-mer size (``k``), ProPhyle version (``prophyle-version``, ``prophyle-revision``, ``prophyle-commit``)
	   * - ``log.txt``
	     - Log
	   * - ``tree.nw``
	     - Phylogenetic tree adjusted for classification
	   * - ``tree.preliminary.nw``
	     - Phylogenetic tree before adjusting


Compressed ProPhyle index for transmission
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A ``.tar.gz`` archive with the following subset of the index files:

	.. list-table:: Compressed ProPhyle index
	   :widths: 5 20
	   :header-rows: 1

	   * - File name
	     - Description
	   * - ``index.fa.amb``
	     - Identical
	   * - ``index.fa.ann``
	     - Identical
	   * - ``index.fa.bwt``
	     - Burrows-Wheeler Transform *without* the OCC table (BWA format, before ``bwa bwtupdate``)
	   * - ``index.fa.kmers.tsv``
	     - Identical
	   * - ``index.json``
	     - Identical
	   * - ``tree.nw``
	     - Identical
	   * - ``tree.preliminary.nw``
	     - Identical
