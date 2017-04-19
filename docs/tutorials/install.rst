.. _install:

Installation
============


Prerequisities
--------------

* GCC 4.8+
* ZLib
* Python 3 with ete3 library
* SamTools


Using Conda (recommended)
-------------------------

Environment installation::

	conda create -c bioconda -n prophyle prophyle

Environment activation::

	source activate prophyle


Using PIP
---------

From PyPI::

	pip install --upgrade prophyle

From Git::

	pip install --upgrade git+https://github.com/karel-brinda/prophyle

From PyPI to the current directory::

	pip install --user prophyle
	export PYTHONUSERBASE=`pwd`
	export PATH=$PATH:`pwd`/bin


Quick test
----------

(small k, subsampled bacterial database)::

	prophyle download bacteria
	prophyle index -k 10 ~/prophyle/test_bacteria.nw test_idx
	prophyle classify test_idx reads.fq > result.sam
