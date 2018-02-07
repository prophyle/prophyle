.. _requirements:

Requirements
============

.. glossary::

	Operating system
		ProPhyle requires a Unix-like system with Python 3.3+ and GCC 4.8+.

	Memory
		Memory consumption can be roughly estimated as the total length of all
		contigs obtained from propagation multiplied by two. Such length approximately
		corresponds to the size of the index.fa file. For instance, the RefSeq bacterial
		database requires 14.2 GB.

	CPU
		ProPhyle does not have any particular requirements on the number of CPU cores.

	Disk space
		The required disk space can be roughly estimated as the size of the original FASTA
		files multiplied by five.
