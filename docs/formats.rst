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

NHX trees
^^^^^^^^^

`New Hampshire X Format <https://sites.google.com/site/cmzmasek/home/software/forester/nhx>`_
is parsed using the `ETE3 library <http://etetoolkit.org/>`_  (see specification of `format 1 <http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-and-writing-newick-trees>`_).

Sequences
^^^^^^^^^

Input sequences can be provided in the FASTA or FASTQ formats. Any non-ACGT characters are treated as
unknown nucleotides. All k-mers containing an unknown nucleotide are discarded.


Assignments
-----------

Read assignments in SAM/BAM
^^^^^^^^^^^^^^^^^^^^^^^^^^^


`SAM format <http://samtools.github.io/hts-specs/>`_ 


Read assignments in a Kraken-like format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Kraken <https://ccb.jhu.edu/software/kraken/MANUAL.html#output-format>`_


Abundances
----------

Abundances in a Kraken report
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`kraken-report <https://ccb.jhu.edu/software/kraken/MANUAL.html#sample-reports>`_ format:

	1. Percentage of reads covered by the clade rooted at this taxon
	2. Number of reads covered by the clade rooted at this taxon
	3. Number of reads assigned directly to this taxon
	4. A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are simply '-'.
	5. NCBI taxonomy ID
	6. indented scientific name


Abundances in a Metaphlan2 report
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`MetaPhlAn2 <https://bitbucket.org/biobakery/biobakery/wiki/metaphlan2#rst-header-output-files>`_ format:
	1. clades, ranging from taxonomic kingdoms (Bacteria, Archaea, etc.) through species.
	2. The taxonomic level of each clade is prefixed to indicate its level: Kingdom: ``k__``, Phylum: ``p__``, Class: ``c__``, Order: ``o__``, Family: ``f__``, Genus: ``g__``, Species: ``s__``.

Since sequence-based profiling is relative and does not provide absolute cellular abundance measures,
	clades are hierarchically summed. Each level will sum to 100%; that is, the sum of all kindom-level
	clades is 100%, the sum of all genus-level clades (including unclassified) is also 100%, and so forth.
	OTU equivalents can be extracted by using only the species-level ``s__`` clades from this file
	(again, making sure to include clades unclassified at this level).


Abundances in a Centrifuge report
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Centrifuge <https://ccb.jhu.edu/software/centrifuge/manual.shtml#centrifuge-summary-output-the-default-filename-is-centrifuge_report.tsv>`_ format::

	#name                                                           taxID   taxRank    kmerCount   numReads   numUniqueReads   abundance
	Wigglesworthia glossinidia endosymbiont of Glossina brevipalpis 36870   leaf       703004      5981.37    5964             0

1. name of a genome, or the name corresponding to a taxonomic ID (the second column) at a rank higher than the strain (e.g., Wigglesworthia glossinidia endosymbiont of Glossina brevipalpis).
2. taxonomic ID (e.g., 36870).
3. taxonomic rank (e.g., leaf).
4. number of k-mers propagated up to the node (e.g., 703004).
5. number of reads classified to this node including multi-classified reads (divided by the number of assignments, e.g., 5981.37).
6. number of reads uniquely classified to this genomic sequence (e.g., 5964).
7. not used yet.



Internal ProPhyle formats
-------------------------

ProPhyle Index
^^^^^^^^^^^^^^

Directory with the following files:

* ``index.fa`` – assembled contigs, name of sequences are of the following format: ``[node_name]@c[contig_number]``
* ``index.fa.amb`` – list of ambiguous nucleotides, no values
* ``index.fa.ann`` – list of contigs and their starting positions in the master string
* ``index.fa.[k].klcp`` – k-LCP array
* ``index.fa.bwt`` – Burrows-Wheeler Transform of the master string (merged sequences + reverse completement) + OCC table (BWA format)
* ``index.fa.kmers.tsv`` – k-mer statistics, format: ``[node_name].[full|reduced].fa	[#kmers]``, where ``full`` refers to all associated k-mers and ``reduced`` to represented k-mers 
* ``index.fa.pac`` – packed sequences (BWA format)
* ``index.fa.sa`` – sampled suffix array (BWA format)
* ``index.json`` – JSON dictionary with a k-mer size (``k``) and ProPhyle version specification (``prophyle-version``, ``prophyle-revision``, ``prophyle-commit``)
* ``log.txt``
* ``tree.nw`` – phylogenetic tree adjusted for classification
* ``tree.preliminary.nw`` – phylogenetic tree before adjusting


Compressed index for transmission
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Archive ``[name].tar.gz`` with the following subset of index files:

* ``index.fa.amb``
* ``index.fa.ann``
* ``index.fa.bwt`` – Burrows-Wheeler Transform *without* the OCC table
* ``index.fa.kmers.tsv``
* ``index.json``
* ``tree.nw``
* ``tree.preliminary.nw``
