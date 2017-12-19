.. _install:

Installing ProPhyle
===================


ProPhyle is written in Python, C and C++. It is distributed as a Python package
and all C/C++ programs are compiled upon the first execution
of the main program. ProPhyle requires a Unix operating system and the following dependencies:

* Python 3.3+ with the `ETE3 library <http://etetoolkit.org/>`_
* GCC 4.8+
* ZLib

There are multiple ways of installation:

.. contents::
	:depth: 1
	:local:
	:backlinks: none


Installing ProPhyle using Bioconda
----------------------------------

To set-up Bioconda, install
`Miniconda <https://conda.io/miniconda.html>`_
or another Conda distribution, and
add `Bioconda channels <https://bioconda.github.io/>`_::

	$ conda config --add channels defaults
	$ conda config --add channels conda-forge
	$ conda config --add channels bioconda

You may either create a separate Conda environment
(which is the recommended approach)::

	$ conda create -n prophyle prophyle
	$ source activate prophyle

or install ProPhyle directly to your main environment::

	$ conda install prophyle


Installing ProPhyle using PIP
-----------------------------

There three options of installing ProPhyle using PIP.

1) From PyPI::

	$ pip install -U prophyle

2) From Git::

	$ pip install -U git+https://github.com/prophyle/prophyle

3) From PyPI to the current directory::

	$ pip install --user prophyle
	$ export PYTHONUSERBASE=`pwd`
	$ export PATH=$PATH:`pwd`/bin



Running ProPhyle directly from the repository
---------------------------------------------

It also is possible to run ProPhyle directly from the repository, by calling
the main script::

	$ prophyle/prophyle/prophyle.py

or its alias::

	$ prophyle/prophyle/prophyle

ProPhyle will then automatically adjust all paths of the auxiliary programs.

Note that ProPhyle uses submodules, therefore the repository needs to
be clonned with the `--recursive` option::

    $ git clone --recursive http://github.com/prophyle/prophyle


Adjusting path
~~~~~~~~~~~~~~

The ProPhyle path can be prepended to the `$PATH` variable so that ProPhyle
can be executed in the same way as if it was installed using PIP::

    $ export PATH=$(pwd)/prophyle/prophyle:$PATH


Installing dependencies
~~~~~~~~~~~~~~~~~~~~~~~

When run from the repository,
some of the ProPhyle dependencies, listed in `requirements.txt`, might
be missing in the system.
It is possible to install them either using BioConda::

    $ cat prophyle/requirements.txt | perl -pe 's/==.*//g' | xargs conda install

or using PIP::

    $ cat prophyle/requirements.txt | xargs pip install
