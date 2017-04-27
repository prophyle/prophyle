.. _install:

Installation
============

Prerequisities
--------------

ProPhyle is written in Python and C++. It is distributed as a Python package
and all C++ auxiliary programs are compiled upon the first execution of the main program.
ProPhyle requires the following dependencies:

* Python 3 with ETE3
* GCC 4.8+
* ZLib


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

To quickly test if ProPhyle has been installed correctly, you can
create a small index with a small k-mer length::

	prophyle download bacteria
	prophyle index -k 10 ~/prophyle/test_bacteria.nw test_idx
