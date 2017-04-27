import setuptools

import os
import sys

if sys.version_info < (3,2):
	sys.exit('Minimum supported Python version is 3.2')

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
	long_description = f.read()


# Get the current version
exec(open("prophyle/version.py").read())

setuptools.setup(
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
		'Programming Language :: Python :: 3 :: Only',
		'Operating System :: Unix',
		'Environment :: Console',
		'License :: OSI Approved :: MIT License',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Bio-Informatics',
	],

	keywords='metagenomics classification NGS',

	packages = ["prophyle"],

	install_requires=['ete3', 'wheel', 'bitarray', 'psutil', 'argparse'],

	package_data={
		'prophyle': [
			'Makefile',
			'*.py',
			'prophyle_assembler/*.cpp',
			'prophyle_assembler/*.h',
			'prophyle_assembler/Makefile',
			'prophyle_index/*.c',
			'prophyle_index/*.h',
			'prophyle_index/Makefile',
			'prophyle_index/bwa/*.c',
			'prophyle_index/bwa/*.h',
			'prophyle_index/bwa/Makefile',
			'trees/*.nw',
		],
	},

	entry_points={
			'console_scripts': [
				'prophyle = prophyle.prophyle:main',
				'prophyle_propagation_makefile.py = prophyle.prophyle_propagation_makefile:main',
				'prophyle_merge_fa.py = prophyle.prophyle_merge_fa:main',
				'prophyle_test_tree.py = prophyle.prophyle_test_tree:main',
				'prophyle_assignment.py = prophyle.prophyle_assignment:main',
				'prophyle_merge_trees.py = prophyle.prophyle_merge_trees:main',
				'prophyle_ncbi_tree.py = prophyle.prophyle_ncbi_tree:main',
			]
	},
)
