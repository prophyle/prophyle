#import ez_setup
#ez_setup.use_setuptools()

from setuptools import setup
from setuptools import Extension
from setuptools import find_packages

import os
import subprocess


from distutils.command.build_ext import build_ext as _build_ext

prophyle_assembler_mod = Extension(
    "fake_extension",
    ["src/fake_extension.c"],
)

class build_ext(_build_ext):
	def run(self):
		subprocess.call('make -C src', shell=True)
		_build_ext.run(self)

setup(
    ext_modules=[prophyle_assembler_mod],
	name='prophyle',
	version='0.0.1',
	description='ProPhyle metagenomic classifier',
	packages = find_packages(),
	url='https://github.com/karel-brinda/prophyle',
	license='MIT',
	cmdclass={'build_ext': build_ext},
	classifiers=[
		'Development Status :: 4 - Beta',
		'Topic :: Scientific/Engineering :: Bio-Informatics',
		'Programming Language :: Python :: 3 :: Only',
		'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
	],
)
