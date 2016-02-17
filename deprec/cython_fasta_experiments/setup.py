from distutils.core import setup, Extension
from Cython.Build import cythonize

import os

os.environ["CC"] = "g++"

# First create an Extension object with the appropriate name and sources
ext = Extension(
		name="wrap_fib",
		sources=["cfib.cpp", "wrap_fib.pyx"],
		#extra_compile_args=["-stdlib=libc++"],
		language="c++",
	)

# Use cythonize on the extension object.
setup(ext_modules=cythonize(
		ext,
	),

)

