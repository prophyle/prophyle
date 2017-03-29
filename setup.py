import ez_setup
ez_setup.use_setuptools()

from setuptools import setup
from setuptools import Extension
from subprocess import call


# Only run lib setup when needed, not on every invocation
from distutils.command.build_ext import build_ext as _build_ext


class build_ext(_build_ext):
	"""Specialized Python extension builder."""

	def run(self):
		call('make -C src', shell=True)
		_build_ext.run(self)


#setup(cmdclass={'build_ext': build_ext}, **setup_metadata)


setup(
	name='prophyle',
	version='0.0.1',
	description='ProPhyle metagenomic classifier',
	#packages = find_packages(),
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
