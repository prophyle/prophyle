:tocdepth: 1

.. _install:

Installation and running
========================

Prerequisities
--------------

ProPhyle is written in Python and C++. It is distributed as a Python package
and all C++ auxiliary programs are compiled upon the first execution of the main program.
ProPhyle requires the following dependencies:

* Python 3.3+ with ETE3
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

There three options how to install ProPhyle using PIP.

From PyPI::

	pip install --upgrade prophyle

From Git::

	pip install --upgrade git+https://github.com/karel-brinda/prophyle

From PyPI to the current directory::

	pip install --user prophyle
	export PYTHONUSERBASE=`pwd`
	export PATH=$PATH:`pwd`/bin


Running ProPhyle directly from the Github repo
----------------------------------------------

It also is possible to run ProPhyle directly from the repository, by calling
either the script `prophyle/prophyle/prophyle.py` or its alias `propphyle/prophyle/prophyle`.
ProPhyle will then automatically adjust the paths of the auxiliary
programs (e.g., `prophyle_index`).

ProPhyle uses submodule, therefore the repository should be clonned with the
`--recursive` option, e.g.::

        git clone --recursive http://github.com/karel-brinda/prophyle

The ProPhyle path can be then added to the `$PATH` variable so that ProPhyle
can be executed in the same way as if it was installed using PIP::

        export PATH=$(pwd)/prophyle/prophyle:$PATH

Note that some of the ProPhyle dependencies, listed in `requirements.txt` might
be missing in the system.  It is possible to install them either using BioConda::

        cat requirements.txt | xargs conda install

or using PIP::

        cat requirements.txt | xargs pip install


Quick test
----------

To quickly test if ProPhyle has been installed correctly, you can
create a small index with a small k-mer length::

	prophyle download bacteria
	prophyle index -k 10 -s 0.1 ~/prophyle/bacteria.nw test_idx
