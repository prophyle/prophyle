#! /usr/bin/env python3

import os
import sys

vfn = os.path.join(os.path.dirname(sys.argv[0]), "version.py")

exec(open(vfn).read())

numbers = VERSION.split(".")
numbers[-1] = str(int(numbers[-1]) + 1)

version = ".".join(numbers)

with open(vfn, "w") as f:
    f.write('try:\n')
    f.write('    from __commit import *\n')
    f.write('except ImportError:\n')
    f.write('    pass\n')
    f.write('VERSION = "{}"'.format(version))
