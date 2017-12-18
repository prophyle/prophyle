:tocdepth: 1

ProPhyle
========

ProPhyle is a phylogeny-based metagenomic classifier,
whose indexing strategy is based on the bottom-up propagation of genomes’
k-mers in the tree, assembling contigs at each node and matching using a
full-text search.

Compared to other state-of-the-arts classifiers, ProPhyle provides several unique features:


* **Low memory requirements.** Compared to `Kraken <https://ccb.jhu.edu/software/kraken/>`_, ProPhyle has 9x smaller memory footprint for index construction and 5x smaller footprint for querying.
* **Flexibility.** ProPhyle is easy to use with any user-provided phylogenetic trees and reference genomes.
* **Standard bioinformatics formats.** `Newick/NHX <https://sites.google.com/site/cmzmasek/home/software/forester/nhx>`_ is used for representing phylogenetic trees and `SAM <https://samtools.github.io/hts-specs/SAMv1.pdf>`_ for reporting the assignments.
* **Lossless k-mer indexing.** ProPhyle stores a list of *all* genomes containing a *k*-mer.  It can be, therefore, accurate even with trees containing similar genomes (e.g, phylogenetic trees for a single species).
* **Deterministic behavior.** ProPhyle is a fully deterministic classifier with a mathematically well-defined behavior.



Issues
------

Please use `Github issues <https://github.com/prophyle/prophyle/issues>`_.


Changelog
---------

See `Releases <https://github.com/prophyle/prophyle/releases>`_.


Licence
-------

`MIT <https://github.com/prophyle/prophyle/blob/master/LICENSE.txt>`_


Authors
-------

`Karel Brinda <http://brinda.cz>`_ <kbrinda@hsph.harvard.edu>

Kamil Salikhov <kamil.salikhov@univ-mlv.fr>

Simone Pignotti <pignottisimone@gmail.com>

`Gregory Kucherov <http://www-igm.univ-mlv.fr/~koutcher/>`_ <gregory.kucherov@univ-mlv.fr>




.. _`GitHub repository`: http://github.com/karel-brinda/ProPhyle
.. _`Bug reporting`: http://github.com/karel-brinda/ProPhyle/issues
.. _`Contact`: 	kbrinda@hsph.harvard.edu


.. toctree::
	:hidden:
	:titlesonly:
	:name: mastertoc

	tutorial
	example
	requirements
	install
	databases_predefined
	databases_custom
	classification
	formats
	reference
