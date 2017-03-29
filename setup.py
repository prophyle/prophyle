import os
import shutil
import sys
import glob
from setuptools import setup, find_packages, Extension

prophyle_assembler_c_files = ["src/prophyle-assembler/prophyle-assembler.cpp"]

prophyle_assembler_mod = Extension(
    "prophyle-assembler",
    prophyle_assembler_c_files,
	extra_compile_args=['-std=c++11', '-v', '-mmacosx-version-min=10.9'],
    extra_link_args=['-lz'],
    include_dirs=[os.path.join('src', 'prophyle-assembler')],
	language = "c++",
)

setup(
    ext_modules=[prophyle_assembler_mod],
    name='prophyle',
    version='0.0.1',
    description='ProPhyle metagenomic classifier',
    packages = find_packages(),
    url='https://github.com/karel-brinda/prophyle',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
)
