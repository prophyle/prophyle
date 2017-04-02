from setuptools import setup
from setuptools import Extension
from setuptools import find_packages

from distutils.command.build_ext import build_ext as _build_ext

import os
import subprocess
import sys

if sys.version_info < (3,2):
	sys.exit('Minimum supported Python version is 3.2')

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
	long_description = f.read()


# Get version
exec(open("prophyle/version.py").read())


## Running make deprec
#fake_ext = Extension(
#	"prophyle.fake_extension",
#	["prophyle/fake_extension.c"],
#)

class build_ext(_build_ext):
	def run(self):
		##
		## Running make deprecated, reasons:
		## - does not work with PIP
		## - when executables present on list of package files, they get uploaded to PyPI
		##
		#subprocess.call('make -C prophyle', shell=True)
		_build_ext.run(self)

setup(
	name='prophyle',

	version=VERSION,

	description='ProPhyle metagenomic classifier',

	long_description=long_description,

	url='https://github.com/karel-brinda/prophyle',

	author='Karel Brinda, Kamil Salikhov, Simone Pignotti, Gregory Kucherov',
	author_email='kbrinda@hsph.harvard.edu, salikhov.kamil@gmail.com, pignottisimone@gmail.com, gregory.kucherov@univ-mlv.fr',

	license='MIT',

	classifiers=[
		'Development Status :: 4 - Beta',
		'Topic :: Scientific/Engineering :: Bio-Informatics',
		'Programming Language :: Python :: 3 :: Only',
		'License :: OSI Approved :: MIT License',
	],

	keywords='metagenomics classification NGS',

	packages = find_packages(),

	install_requires=['ete3', 'numpy', 'wheel', 'six', 'scipy', 'bitarray'],

	package_data={
		'prophyle': [
			'Makefile',
			'*.py',
			'prophyle-assembler/*.cpp',
			'prophyle-assembler/*.h',
			'prophyle-assembler/Makefile',
			'prophyle-index/*.c',
			'prophyle-index/*.h',
			'prophyle-index/Makefile',
			'prophyle-index/bwa/*.c',
			'prophyle-index/bwa/*.h',
			'prophyle-index/bwa/Makefile',
			'trees/*.nw',
		],
	},

	## Running make deprec
	#ext_modules=[fake_ext],

	entry_points={
			'console_scripts': [
			'prophyle = prophyle.prophyle:main',
			'newick2makefile.py = prophyle.newick2makefile:main',
			'create_final_fasta.py = prophyle.create_final_fasta:main',
			'assignment.py = prophyle.assignment:main',
		]
	},

	cmdclass={'build_ext': build_ext},
)
