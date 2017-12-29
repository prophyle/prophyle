.. _requirements:

Requirements
============

.. glossary::

	Operating system
		ProPhyle requires a Unix-like system with Python 3.3+ and GCC 4.8+.

	Memory
		Memory consumption is slightly lower for index construction than classification
		and depends on the source data. It can be estimated as ....
		For an existing index, the memory footprint can be estimated as the sum
		of .... If k-LCP array is not used, the footprint will be smaller by its size.

	CPU
		ProPhyle does not have any particular requirements on the number of CPU cores.
		In the current version, it cannot efficiently use multiple cores.

	Disk space
		ProPhyle requires approximately ... for index construction.
		The resulting index usually occupies approximately ....

	Time
		...
